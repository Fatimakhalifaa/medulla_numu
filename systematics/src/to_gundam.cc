/**
 * @file to_gundam.cc
 * @brief Convert a (Medulla + run_systematics) output file into GUNDAM input.
 * @details One pass per [[tree]]:
 *   - copies every branch of the event tree by type (int / float / double /
 *     Long64 / bool), and writes the GUNDAM flags as ints
 *     (Run, Subrun, Evt, cut_type, is_nu, is_data, category; same rules as
 *     the previous to_gundam);
 *   - merges the systematics trees <origin>_<type>Tree into the same tree,
 *     entry by entry (Run/Subrun/Evt are checked on every entry):
 *       multisigma, variation, NuMIfluxsim -> one TClonesArray(TGraph) per dial
 *       multisim                           -> vector branches (or a clone)
 *     vector<float> and vector<double> are both accepted;
 *   - adds the phase-space fill branches (nu_mu only, from the notebook
 *     7.17 table): fill_weight (weight option), phase_space_pass (cut
 *     option), fill_window_pass, fill_topology, fill_frac, fill_u,
 *     fill_elep, fill_evis. Without a table they are written as 1 (nu_e);
 *   - copies POT / Livetime from the origin directory, optionally scaled
 *     (failed-SPINE-spill correction).
 * Usage: to_gundam <config.toml> [input.root] [output.root]
 */
#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <TBranch.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>

#include "configuration.h"
#include "gundam_fill.h"

// ---------------------------------------------------------------------------
// Config helpers (optional keys fall back to a default)
// ---------------------------------------------------------------------------
static std::string opt_string(cfg::ConfigurationTable & t, const std::string & key, const std::string & def)
{
    try { std::string v = t.get_string_field(key); return v.empty() ? def : v; }
    catch(...) { return def; }
}
static bool opt_bool(cfg::ConfigurationTable & t, const std::string & key, bool def)
{
    try { return t.get_bool_field(key); }
    catch(...) { return def; }
}
static double opt_double(cfg::ConfigurationTable & t, const std::string & key, double def)
{
    const std::string s = opt_string(t, key, "");
    if(s.empty()) return def;
    return std::stod(s);
}

// ---------------------------------------------------------------------------
// Directories
// ---------------------------------------------------------------------------
static std::string trim_slashes(std::string p)
{
    while(!p.empty() && p.back() == '/') p.pop_back();
    while(!p.empty() && p.front() == '/') p.erase(p.begin());
    return p;
}
static std::string dirname_of(const std::string & path)
{
    const std::string p = trim_slashes(path);
    const size_t i = p.rfind('/');
    return i == std::string::npos ? std::string() : p.substr(0, i);
}
static TDirectory * make_dirs(TDirectory * top, const std::string & path)
{
    TDirectory * d = top;
    std::string p = trim_slashes(path);
    size_t start = 0;
    while(start < p.size())
    {
        size_t end = p.find('/', start);
        if(end == std::string::npos) end = p.size();
        const std::string part = p.substr(start, end - start);
        if(!part.empty())
        {
            TDirectory * next = d->GetDirectory(part.c_str());
            if(!next) next = d->mkdir(part.c_str());
            if(!next) throw std::runtime_error("cannot create directory " + path);
            d = next;
        }
        start = end + 1;
    }
    return d;
}

// ---------------------------------------------------------------------------
// Scalar branch buffer (input type kept on output)
// ---------------------------------------------------------------------------
struct ScalarBuffer
{
    std::string name;
    char type = 'D';        // I, F, D, L (Long64), O (bool)
    Int_t i = 0; Float_t f = 0; Double_t d = 0; Long64_t l = 0; Bool_t o = false;

    bool bind_input(TTree * tree, TBranch * br)
    {
        name = br->GetName();
        TObjArray * leaves = br->GetListOfLeaves();
        if(!leaves || leaves->GetEntries() != 1) return false;
        TLeaf * leaf = (TLeaf *)leaves->At(0);
        if(leaf->GetLen() != 1) return false;   // arrays are not scalars
        const std::string t = leaf->GetTypeName();
        void * addr = nullptr;
        if(t == "Int_t" || t == "int") { type = 'I'; addr = &i; }
        else if(t == "Float_t" || t == "float") { type = 'F'; addr = &f; }
        else if(t == "Double_t" || t == "double") { type = 'D'; addr = &d; }
        else if(t == "Long64_t" || t == "Long_t" || t == "long") { type = 'L'; addr = &l; }
        else if(t == "Bool_t" || t == "bool") { type = 'O'; addr = &o; }
        else return false;
        tree->SetBranchAddress(name.c_str(), addr);
        return true;
    }
    void bind_output(TTree * tree, const std::string & out_name = "")
    {
        const std::string n = out_name.empty() ? name : out_name;
        void * addr = (type == 'I') ? (void *)&i : (type == 'F') ? (void *)&f : (type == 'L') ? (void *)&l
                    : (type == 'O') ? (void *)&o : (void *)&d;
        tree->Branch(n.c_str(), addr, (n + "/" + type).c_str());
    }
    double value() const
    {
        switch(type)
        {
            case 'I': return i;
            case 'F': return f;
            case 'L': return double(l);
            case 'O': return o ? 1.0 : 0.0;
            default: return d;
        }
    }
};

// ---------------------------------------------------------------------------
// vector<float> / vector<double> input
// ---------------------------------------------------------------------------
struct VecIn
{
    std::string name;
    bool is_float = false;
    std::vector<float> * vf = nullptr;
    std::vector<double> * vd = nullptr;

    bool bind(TTree * tree, const std::string & n)
    {
        name = n;
        TBranch * br = tree->GetBranch(n.c_str());
        if(!br) return false;
        const std::string cls = br->GetClassName();
        // objects are allocated here, so the input and output branches share them
        if(cls == "vector<float>") { is_float = true; vf = new std::vector<float>(); tree->SetBranchAddress(n.c_str(), &vf); return true; }
        if(cls == "vector<double>") { is_float = false; vd = new std::vector<double>(); tree->SetBranchAddress(n.c_str(), &vd); return true; }
        std::cerr << "  [warning] " << n << " has type '" << cls << "', skipped\n";
        return false;
    }
    size_t size() const { return is_float ? (vf ? vf->size() : 0) : (vd ? vd->size() : 0); }
    double at(size_t k) const { return is_float ? double((*vf)[k]) : (*vd)[k]; }
};

// ---------------------------------------------------------------------------
// One systematics tree merged into the output tree
// ---------------------------------------------------------------------------
struct SystTree
{
    std::string type;                // multisigma, variation, NuMIfluxsim, multisim
    TTree * tree = nullptr;
    Int_t run = 0, subrun = 0, evt = 0;
    bool has_id = false;
    // TGraph dials
    std::deque<VecIn> weights, nsigmas;
    std::vector<TClonesArray *> graphs;
    // multisim vectors (copied as they are)
    std::deque<VecIn> vectors;
    long n_nonfinite = 0;
};

static const std::vector<std::string> kMorphDials = {
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad"};

static bool bind_id(TTree * t, Int_t & run, Int_t & subrun, Int_t & evt)
{
    for(const char * b : {"Run", "Subrun", "Evt"})
    {
        TBranch * br = t->GetBranch(b);
        if(!br) return false;
        TLeaf * leaf = (TLeaf *)br->GetListOfLeaves()->At(0);
        const std::string tn = leaf ? leaf->GetTypeName() : "";
        if(tn != "Int_t" && tn != "int") return false;
    }
    t->SetBranchAddress("Run", &run);
    t->SetBranchAddress("Subrun", &subrun);
    t->SetBranchAddress("Evt", &evt);
    return true;
}

/** Set the event's graph: sorted (x, y) points; non-finite weights -> 1. */
static void set_graph(TClonesArray * arr, std::vector<double> & x, std::vector<double> & y, long & n_nonfinite)
{
    std::vector<size_t> order(x.size());
    for(size_t k = 0; k < order.size(); ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
    TGraph * g = (TGraph *)arr->ConstructedAt(0);
    g->Set(int(x.size()));
    for(size_t k = 0; k < order.size(); ++k)
    {
        double w = y[order[k]];
        if(!gfill::finite(w)) { w = 1.0; ++n_nonfinite; }
        g->SetPoint(int(k), x[order[k]], w);
    }
}

// ---------------------------------------------------------------------------
// Fill configuration (branch names in the event tree)
// ---------------------------------------------------------------------------
struct FillConfig
{
    bool enabled = false;
    gfill::FillTable table;
    std::string containment, elep, evis, cos, ehad, rest, emu, vx, vy, vz;
    double emu_scale = 1.0;
};


static int process_tree(cfg::ConfigurationTable & config, cfg::ConfigurationTable & table,
                        TFile * input, TFile * output, const FillConfig & fill, const std::string & multisim_mode,
                        double ms_cap)
{
    const std::string origin = trim_slashes(table.get_string_field("origin"));
    const std::string destination = table.get_string_field("destination");
    const std::string name = table.get_string_field("name");
    const bool tree_is_nu = table.get_bool_field("is_nu");
    const bool tree_is_data = table.get_bool_field("is_data");
    const bool store_syst = opt_bool(table, "gundam_store_syst", false);

    const bool do_fill = fill.enabled && opt_bool(table, "fill", true);   // per-tree switch; off -> pass = weight = 1

    TTree * in = (TTree *)input->Get(origin.c_str());
    if(!in) { std::cerr << "Warning: tree " << origin << " not found, skipped\n"; return 0; }
    std::cout << "Processing " << origin << " (" << in->GetEntries() << " entries) -> " << destination << name << "\n";
    if(fill.enabled && !do_fill) std::cout << "  fill disabled for this tree (phase_space_pass = fill_weight = 1)\n";

    TDirectory * dir = make_dirs(output, destination);

    // ---- POT / Livetime from the origin directory, optionally scaled ----
    const std::string src_dir = dirname_of(origin);
    const double pot_scale = opt_double(table, "pot_scale", 1.0);
    const double live_scale = opt_double(table, "livetime_scale", pot_scale);
    for(const auto & [hname, scale] : {std::make_pair(std::string("POT"), pot_scale), std::make_pair(std::string("Livetime"), live_scale)})
    {
        TH1 * h = (TH1 *)input->Get((src_dir.empty() ? hname : src_dir + "/" + hname).c_str());
        if(!h) { std::cerr << "  [warning] no " << hname << " in " << src_dir << "\n"; continue; }
        TH1 * c = (TH1 *)h->Clone(hname.c_str());
        c->SetDirectory(nullptr);
        if(scale != 1.0) c->Scale(scale);
        dir->WriteTObject(c, hname.c_str(), "Overwrite");
        std::cout << "  " << hname << " = " << c->Integral() << (scale != 1.0 ? " (scaled by " + std::to_string(scale) + ")" : std::string()) << "\n";
        delete c;
    }

    dir->cd();
    TTree * out = new TTree(name.c_str(), name.c_str());
    out->SetDirectory(dir);

    // ---- event-tree branches ----
    ScalarBuffer run_b, subrun_b, evt_b, cut_b, nu_b, data_b, cat_b;
    bool has_run = false, has_subrun = false, has_evt = false, has_cut = false, has_nu = false, has_data = false, has_cat = false;
    std::deque<ScalarBuffer> others;
    std::map<std::string, const ScalarBuffer *> by_name;
    TObjArray * branches = in->GetListOfBranches();
    for(int k = 0; k < branches->GetEntries(); ++k)
    {
        TBranch * br = (TBranch *)branches->At(k);
        const std::string b = br->GetName();
        ScalarBuffer * target = nullptr;
        if(b == "Run" || b == "run") { target = &run_b; has_run = true; }
        else if(b == "Subrun" || b == "subrun") { target = &subrun_b; has_subrun = true; }
        else if(b == "Evt" || b == "evt" || b == "Event" || b == "event") { target = &evt_b; has_evt = true; }
        else if(b == "reco_cut_type" || b == "cut_type") { target = &cut_b; has_cut = true; }
        else if(b == "reco_is_nu" || b == "is_nu") { target = &nu_b; has_nu = true; }
        else if(b == "reco_is_data" || b == "is_data") { target = &data_b; has_data = true; }
        else if(b == "true_category" || b == "category") { target = &cat_b; has_cat = true; }
        if(target)
        {
            if(!target->bind_input(in, br)) throw std::runtime_error("unsupported type for " + b);
            continue;
        }
        others.emplace_back();
        if(!others.back().bind_input(in, br))
        {
            std::cerr << "  [warning] branch " << b << " is not a supported scalar, skipped\n";
            others.pop_back();
            in->SetBranchStatus(b.c_str(), false);
            continue;
        }
        others.back().bind_output(out);
        by_name[b] = &others.back();
    }

    Int_t run = 0, subrun = 0, evt = 0, cut_type = 0, is_nu = 0, is_data = 0, category = 0;
    out->Branch("Run", &run, "Run/I");
    out->Branch("Subrun", &subrun, "Subrun/I");
    out->Branch("Evt", &evt, "Evt/I");
    out->Branch("cut_type", &cut_type, "cut_type/I");
    out->Branch("is_nu", &is_nu, "is_nu/I");
    out->Branch("is_data", &is_data, "is_data/I");
    out->Branch("category", &category, "category/I");

    // ---- fill branches ----
    Int_t f_pass = 1, f_win = 1, f_topo = -1;
    Double_t f_weight = 1, f_frac = 1, f_u = -1, f_elep = NAN, f_evis = NAN;
    out->Branch("phase_space_pass", &f_pass, "phase_space_pass/I");
    out->Branch("fill_weight", &f_weight, "fill_weight/D");
    out->Branch("fill_window_pass", &f_win, "fill_window_pass/I");
    out->Branch("fill_topology", &f_topo, "fill_topology/I");
    out->Branch("fill_frac", &f_frac, "fill_frac/D");
    out->Branch("fill_u", &f_u, "fill_u/D");
    out->Branch("fill_elep", &f_elep, "fill_elep/D");
    out->Branch("fill_evis", &f_evis, "fill_evis/D");
    std::vector<const ScalarBuffer *> fin;   // containment, elep, evis, cos, ehad, rest, emu, vx, vy, vz
    if(do_fill)
    {
        for(const std::string & b : {fill.containment, fill.elep, fill.evis, fill.cos, fill.ehad, fill.rest, fill.emu, fill.vx, fill.vy, fill.vz})
        {
            auto it = by_name.find(b);
            if(it != by_name.end()) { fin.push_back(it->second); continue; }
            if(b != fill.emu) throw std::runtime_error("fill: branch '" + b + "' not found in " + origin);
            std::cerr << "  [warning] fill: no '" << b << "' in " << origin << ": exiting events fail the fill\n";
            fin.push_back(nullptr);
        }
    }

    // ---- systematics trees ----
    std::vector<std::unique_ptr<SystTree>> systs;
    if(store_syst)
    {
        std::map<std::string, std::vector<std::string>> sys_names;
        for(cfg::ConfigurationTable & s : config.get_subtables("sys"))
            sys_names[s.get_string_field("type")].push_back(s.get_string_field("name"));
        const std::string prefix = opt_string(table, "syst_prefix", origin);
        for(const std::string type : {"multisigma", "variation", "NuMIfluxsim", "multisim"})
        {
            if(type == "multisim" && multisim_mode != "merge") continue;
            const std::string path = prefix + "_" + type + "Tree";
            TTree * st = (TTree *)input->Get(path.c_str());
            if(!st) { std::cout << "  no " << path << "\n"; continue; }
            if(st->GetEntries() != in->GetEntries())
                throw std::runtime_error(path + " has " + std::to_string(st->GetEntries()) + " entries, event tree has " + std::to_string(in->GetEntries()));
            auto sp = std::make_unique<SystTree>();
            sp->type = type; sp->tree = st;
            sp->has_id = bind_id(st, sp->run, sp->subrun, sp->evt);
            if(!sp->has_id) std::cerr << "  [warning] " << path << ": no int Run/Subrun/Evt, entry order is assumed\n";
            if(type == "multisim")
            {
                // listed [[sys]] of type multisim, or every vector branch
                std::vector<std::string> names = sys_names["multisim"];
                if(names.empty())
                {
                    TObjArray * bl = st->GetListOfBranches();
                    for(int k = 0; k < bl->GetEntries(); ++k)
                    {
                        const std::string b = ((TBranch *)bl->At(k))->GetName();
                        if(b != "Run" && b != "Subrun" && b != "Evt") names.push_back(b);
                    }
                }
                for(const std::string & n : names)
                {
                    if(by_name.count(n)) { std::cerr << "  [warning] " << n << " already in the event tree, skipped\n"; continue; }
                    sp->vectors.emplace_back();
                    if(!sp->vectors.back().bind(st, n)) { sp->vectors.pop_back(); continue; }
                    VecIn & v = sp->vectors.back();
                    if(v.is_float) out->Branch(n.c_str(), &v.vf);
                    else out->Branch(n.c_str(), &v.vd);
                }
                std::cout << "  " << path << ": " << sp->vectors.size() << " multisim vectors\n";
            }
            else
            {
                for(const std::string & n : sys_names[type])
                {
                    if(!st->GetBranch(n.c_str())) { std::cerr << "  [warning] dial " << n << " not in " << path << "\n"; continue; }
                    if(!st->GetBranch((n + "_sigma").c_str())) { std::cerr << "  [warning] " << n << "_sigma not in " << path << "\n"; continue; }
                    sp->weights.emplace_back();
                    sp->nsigmas.emplace_back();
                    if(!sp->weights.back().bind(st, n) || !sp->nsigmas.back().bind(st, n + "_sigma"))
                    {
                        st->SetBranchStatus(n.c_str(), false);
                        st->SetBranchStatus((n + "_sigma").c_str(), false);
                        sp->weights.pop_back(); sp->nsigmas.pop_back();
                        continue;
                    }
                }
                sp->graphs.resize(sp->weights.size(), nullptr);   // sized once: the branches keep &graphs[k]
                for(size_t k = 0; k < sp->weights.size(); ++k)
                {
                    sp->graphs[k] = new TClonesArray("TGraph", 1);
                    out->Branch(sp->weights[k].name.c_str(), &sp->graphs[k], 32000, -1);
                }
                std::cout << "  " << path << ": " << sp->weights.size() << " dials (of " << sys_names[type].size() << " in [[sys]])\n";
            }
            systs.push_back(std::move(sp));
        }
    }

    // ---- entry loop ----
    long n_pass = 0, n_exit = 0, n_exit_pass = 0, n_frac_gt1 = 0;
    double sum_w = 0;
    std::vector<double> xs, ys;
    const Long64_t n = in->GetEntries();
    for(Long64_t e = 0; e < n; ++e)
    {
        in->GetEntry(e);
        run = has_run ? Int_t(run_b.value()) : 0;
        subrun = has_subrun ? Int_t(subrun_b.value()) : 0;
        evt = has_evt ? Int_t(evt_b.value()) : 0;
        cut_type = has_cut ? Int_t(cut_b.value()) : 0;
        if(tree_is_nu) is_nu = 1;
        else is_nu = has_nu ? Int_t(nu_b.value()) : 0;
        if(tree_is_data) { is_data = 1; category = 10; }
        else
        {
            is_data = has_data ? Int_t(data_b.value()) : 0;
            category = has_cat ? Int_t(cat_b.value()) : 0;
        }

        if(do_fill)
        {
            gfill::FillInput fi;
            fi.containment = fin[0]->value(); fi.elep = fin[1]->value(); fi.evis = fin[2]->value();
            fi.cos = fin[3]->value(); fi.ehad = fin[4]->value(); fi.rest = fin[5]->value();
            fi.emu_exiting = fin[6] ? fin[6]->value() * fill.emu_scale : NAN;
            fi.vx = fin[7]->value(); fi.vy = fin[8]->value(); fi.vz = fin[9]->value();
            const gfill::FillResult r = gfill::evaluate(fill.table, fi);
            f_pass = r.pass; f_win = r.window_pass; f_topo = r.topology; f_frac = r.frac;
            f_weight = r.weight; f_u = r.u; f_elep = r.elep; f_evis = r.evis;
            n_pass += r.pass; sum_w += r.weight;
            if(r.topology == 1) { ++n_exit; n_exit_pass += r.pass; }
            if(r.window_pass && r.frac > 1.0) ++n_frac_gt1;
        }

        for(auto & sp : systs)
        {
            sp->tree->GetEntry(e);
            if(sp->has_id && has_run && has_subrun && has_evt && (sp->run != run || sp->subrun != subrun || sp->evt != evt))
                throw std::runtime_error(sp->type + " tree is out of step with " + origin + " at entry " + std::to_string(e) +
                                         " (" + std::to_string(run) + "/" + std::to_string(subrun) + "/" + std::to_string(evt) + " vs " +
                                         std::to_string(sp->run) + "/" + std::to_string(sp->subrun) + "/" + std::to_string(sp->evt) + ")");

            if(ms_cap > 0)
                for(VecIn & v : sp->vectors)
                {
                    if(v.is_float) { for(float & w : *v.vf) if(!gfill::finite(w) || w < 0 || w > ms_cap) { w = 1.0f; ++sp->n_nonfinite; } }
                    else           { for(double & w : *v.vd) if(!gfill::finite(w) || w < 0 || w > ms_cap) { w = 1.0; ++sp->n_nonfinite; } }
                }            
                                         
            for(size_t k = 0; k < sp->weights.size(); ++k)
            {
                const VecIn & wv = sp->weights[k];
                const VecIn & sv = sp->nsigmas[k];
                xs.clear(); ys.clear();
                const bool morph = sp->type == "multisigma" &&
                                   std::find(kMorphDials.begin(), kMorphDials.end(), wv.name) != kMorphDials.end();
                if(morph)
                {
                    // 5-point symmetric spline (same as the previous to_gundam)
                    double w1 = wv.size() ? wv.at(0) : 1.0;
                    if(!gfill::finite(w1) || !tree_is_nu) w1 = 1.0;
                    const double w05 = 1.0 + 0.5 * (w1 - 1.0);
                    xs = {-1, -0.5, 0, 0.5, 1};
                    ys = {w1, w05, 1, w05, 1};
                }
                else
                {
                    const size_t m = std::min(wv.size(), sv.size());
                    for(size_t q = 0; q < m; ++q) { xs.push_back(sv.at(q)); ys.push_back(tree_is_nu ? wv.at(q) : -5.0); }
                    if(sp->type == "multisigma" || xs.empty())
                    {
                        // weight = 1 at nsigma = 0 (also for events without weights)
                        if(std::find(xs.begin(), xs.end(), 0.0) == xs.end()) { xs.push_back(0.0); ys.push_back(tree_is_nu ? 1.0 : -5.0); }
                    }
                }
                set_graph(sp->graphs[k], xs, ys, sp->n_nonfinite);
            }
        }
        out->Fill();
    }

    for(auto & sp : systs)
        if(sp->n_nonfinite) std::cerr << "  [warning] " << sp->type << ": " << sp->n_nonfinite << " non-finite weights set to 1\n";
    if(do_fill)
    {
        std::cout << "  fill: phase_space_pass " << n_pass << " / " << n << ", sum fill_weight " << sum_w
                  << ", exiting " << n_exit << " (pass " << n_exit_pass << ")\n";
        if(n_frac_gt1) std::cout << "  fill: " << n_frac_gt1 << " in-window events with f > 1 (weight > 1, cut keeps them all)\n";
    }

    dir->cd();
    out->Write("", TObject::kOverwrite);
    std::cout << "  wrote " << out->GetEntries() << " entries\n";
    delete out;
    for(auto & sp : systs) for(TClonesArray * a : sp->graphs) delete a;
    return 1;
}

int main(int argc, char * argv[])
{
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <config_file.toml> [input_file.root] [output_file.root]\n";
        return 1;
    }
    try
    {
        cfg::ConfigurationTable config;
        config.set_config(argv[1]);
        std::vector<cfg::ConfigurationTable> tables = config.get_subtables("tree");

        const std::string input_name = argc > 2 ? std::string(argv[2]) : config.get_string_field("output.path");
        const std::string output_name = argc > 3 ? std::string(argv[3])
                                      : opt_string(config, "gundam.output", "gundam_" + input_name.substr(input_name.rfind('/') + 1));
        if(output_name == input_name) throw std::runtime_error("output file is the input file");
        const std::string multisim_mode = opt_string(config, "gundam.multisim", "merge");   // merge | clone | none
        const double ms_cap = opt_double(config, "gundam.multisim_cap", 0.0);              // 0 = off; weights > cap (or non-finite/negative) -> 1

        FillConfig fill;
        const std::string table_path = opt_string(config, "gundam.fill_table", "");
        if(!table_path.empty())
        {
            fill.enabled = true;
            fill.table = gfill::FillTable::load(table_path);
            fill.containment = opt_string(config, "gundam.fill_containment", "reco_leading_muon_containment");
            fill.elep = opt_string(config, "gundam.fill_elep", "reco_leptonic_visible_energy");
            fill.evis = opt_string(config, "gundam.fill_evis", "reco_visible_energy");
            fill.cos = opt_string(config, "gundam.fill_cos", "reco_leading_muon_cos_theta_mu");
            fill.ehad = opt_string(config, "gundam.fill_ehad", "reco_hadronic_visible_energy");
            fill.rest = opt_string(config, "gundam.fill_rest", "reco_hadronic_visible_energy_rest");
            fill.emu = opt_string(config, "gundam.exiting_energy", "true_muon_energy_synth");
            fill.emu_scale = opt_double(config, "gundam.exiting_energy_scale", 1.0);
            fill.vx = opt_string(config, "gundam.fill_vertex_x", "reco_vertex_x");
            fill.vy = opt_string(config, "gundam.fill_vertex_y", "reco_vertex_y");
            fill.vz = opt_string(config, "gundam.fill_vertex_z", "reco_vertex_z");
            std::cout << "Fill table " << table_path << ": " << fill.table.ncell() << " cells, exiting window ["
                      << fill.table.win_x_lo << ", " << fill.table.win_x_hi << ") GeV, E_lep >= " << fill.table.elep_min
                      << ", exiting E_mu from " << fill.emu << " x " << fill.emu_scale << "\n";
        }
        else std::cout << "No gundam.fill_table: phase_space_pass = fill_weight = 1 for every event\n";

        TFile * input = TFile::Open(input_name.c_str(), "READ");
        if(!input || input->IsZombie()) throw std::runtime_error("cannot open input " + input_name);
        TFile * output = TFile::Open(output_name.c_str(), "RECREATE");
        if(!output || output->IsZombie()) throw std::runtime_error("cannot create output " + output_name);
        std::cout << "Input  " << input_name << "\nOutput " << output_name << "\n";

        std::set<std::string> cloned;
        for(cfg::ConfigurationTable & table : tables)
        {
            process_tree(config, table, input, output, fill, multisim_mode, ms_cap);
            if(multisim_mode == "clone" && opt_bool(table, "gundam_store_syst", false))
            {
                const std::string path = opt_string(table, "syst_prefix", trim_slashes(table.get_string_field("origin"))) + "_multisimTree";
                TTree * ms = (TTree *)input->Get(path.c_str());
                if(ms && !cloned.count(path))
                {
                    TDirectory * d = make_dirs(output, table.get_string_field("destination"));
                    d->cd();
                    TTree * c = ms->CloneTree(-1);
                    c->Write("", TObject::kOverwrite);
                    cloned.insert(path);
                    std::cout << "  cloned " << path << "\n";
                }
            }
        }
        output->Close();
        input->Close();
    }
    catch(const std::exception & ex)
    {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
