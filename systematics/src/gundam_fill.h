/**
 * @file gundam_fill.h
 * @brief Phase-space fill (contained + selected exiting nu_mu) for to_gundam.
 * @details No ROOT dependence. Reads the table written by notebook cell 7.17
 * (exiting_fill_table.json) and evaluates, per event:
 *   fill_weight      = cell keep fraction f (0 outside the window)   -> weight option
 *   phase_space_pass = window && u < f, u = vertex_uniform(reco vertex) -> cut option
 * Same binning, hash and conventions as cuts::phase_space_fill_table4_cut
 * (phase_space_fill.h) and exiting_fill_keep() (notebook 7.17).
 */
#ifndef GUNDAM_FILL_H
#define GUNDAM_FILL_H
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace gfill
{
    /** isfinite that survives -Ofast / -ffast-math (the CMake flags), which may fold std::isfinite to true. */
    inline bool finite(double x)
    {
        uint64_t b; std::memcpy(&b, &x, sizeof b);
        return ((b >> 52) & 0x7FFULL) != 0x7FFULL;
    }

    /** Same as utilities::vertex_uniform in phase_space_fill.h (float32 bit patterns, splitmix64). */
    inline double vertex_uniform(float x, float y, float z)
    {
        const float v[3] = {x, y, z};
        uint32_t b[3];
        for(int k = 0; k < 3; ++k) std::memcpy(&b[k], &v[k], sizeof(float));
        uint64_t h = (uint64_t(b[0]) << 32) | uint64_t(b[1]);
        h ^= uint64_t(b[2]) * 0x9E3779B97F4A7C15ULL;
        h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ULL;
        h ^= h >> 27; h *= 0x94D049BB133111EBULL;
        h ^= h >> 31;
        return double(h >> 11) * (1.0 / 9007199254740992.0);
    }

    /** np.searchsorted(edges, x, side="right") - 1, clipped to the first/last bin. */
    inline size_t table_bin(const std::vector<double> & edges, double x)
    {
        long i = long(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
        return size_t(std::clamp<long>(i, 0L, long(edges.size()) - 2));
    }

    // ------------------------------------------------------------------
    // Minimal JSON reader: numbers, strings, bools/null, arrays, objects.
    // Arrays of numbers (any nesting) are flattened in C order.
    // ------------------------------------------------------------------
    struct JsonValue
    {
        enum Kind { kNull, kNumber, kString, kArray, kObject } kind = kNull;
        double num = 0.0;
        std::string str;
        std::vector<JsonValue> arr;
        std::map<std::string, JsonValue> obj;

        void flatten(std::vector<double> & out) const
        {
            if(kind == kNumber) out.push_back(num);
            else if(kind == kArray) for(const auto & v : arr) v.flatten(out);
            else if(kind == kNull) out.push_back(std::nan(""));
            else throw std::runtime_error("JSON: expected a number or an array of numbers");
        }
        std::vector<double> numbers() const { std::vector<double> o; flatten(o); return o; }
        bool has(const std::string & k) const { return kind == kObject && obj.count(k); }
        const JsonValue & at(const std::string & k) const
        {
            auto it = obj.find(k);
            if(it == obj.end()) throw std::runtime_error("JSON: missing key '" + k + "'");
            return it->second;
        }
    };

    class JsonParser
    {
    public:
        explicit JsonParser(const std::string & s) : s_(s) {}
        JsonValue parse() { JsonValue v = value(); ws(); if(p_ != s_.size()) fail("trailing characters"); return v; }
    private:
        const std::string & s_;
        size_t p_ = 0;
        [[noreturn]] void fail(const std::string & m) { throw std::runtime_error("JSON parse error at " + std::to_string(p_) + ": " + m); }
        void ws() { while(p_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[p_]))) ++p_; }
        bool lit(const char * w) { size_t n = std::strlen(w); if(s_.compare(p_, n, w) == 0) { p_ += n; return true; } return false; }
        std::string string()
        {
            std::string o; ++p_;
            while(p_ < s_.size() && s_[p_] != '"')
            {
                if(s_[p_] == '\\') { ++p_; if(p_ >= s_.size()) fail("bad escape"); if(s_[p_] == 'u') { p_ += 4; o += '?'; } else o += s_[p_]; }
                else o += s_[p_];
                ++p_;
            }
            if(p_ >= s_.size()) fail("unterminated string");
            ++p_;
            return o;
        }
        JsonValue value()
        {
            ws();
            if(p_ >= s_.size()) fail("unexpected end");
            JsonValue v;
            char c = s_[p_];
            if(c == '{')
            {
                v.kind = JsonValue::kObject; ++p_; ws();
                if(s_[p_] == '}') { ++p_; return v; }
                while(true)
                {
                    ws(); if(s_[p_] != '"') fail("expected key");
                    std::string k = string(); ws();
                    if(s_[p_] != ':') fail("expected ':'");
                    ++p_;
                    v.obj[k] = value(); ws();
                    if(s_[p_] == ',') { ++p_; continue; }
                    if(s_[p_] == '}') { ++p_; break; }
                    fail("expected ',' or '}'");
                }
            }
            else if(c == '[')
            {
                v.kind = JsonValue::kArray; ++p_; ws();
                if(s_[p_] == ']') { ++p_; return v; }
                while(true)
                {
                    v.arr.push_back(value()); ws();
                    if(s_[p_] == ',') { ++p_; continue; }
                    if(s_[p_] == ']') { ++p_; break; }
                    fail("expected ',' or ']'");
                }
            }
            else if(c == '"') { v.kind = JsonValue::kString; v.str = string(); }
            else if(lit("true")) { v.kind = JsonValue::kNumber; v.num = 1; }
            else if(lit("false")) { v.kind = JsonValue::kNumber; v.num = 0; }
            else if(lit("null")) { v.kind = JsonValue::kNull; }
            else if(lit("NaN")) { v.kind = JsonValue::kNumber; v.num = std::nan(""); }
            else if(lit("Infinity")) { v.kind = JsonValue::kNumber; v.num = 1e300; }
            else if(lit("-Infinity")) { v.kind = JsonValue::kNumber; v.num = -1e300; }
            else
            {
                size_t used = 0;
                try { v.num = std::stod(s_.substr(p_, 64), &used); } catch(...) { fail("bad number"); }
                v.kind = JsonValue::kNumber; p_ += used;
            }
            return v;
        }
    };

    // ------------------------------------------------------------------
    // Fill table
    // ------------------------------------------------------------------
    struct FillTable
    {
        std::vector<std::vector<double>> edges;   // E_lep, cos, E_had, E_had - KE_p
        std::vector<double> frac_contained, frac_exiting;
        double win_x_lo = -1e300, win_x_hi = 1e300;   // exiting E_vis window (proxy), GeV
        bool has_win_c = false;
        double win_c_lo = -1e300, win_c_hi = 1e300;   // optional contained window (already in Medulla)
        double elep_min = -1e300;                        // exiting proxy E_lep threshold, GeV

        size_t ncell() const { size_t n = 1; for(const auto & e : edges) n *= e.size() - 1; return n; }
        size_t cell(double elep, double cos, double ehad, double rest) const
        {
            const double x[4] = {elep, cos, ehad, rest};
            size_t i = 0;
            for(size_t a = 0; a < 4; ++a) i = i * (edges[a].size() - 1) + table_bin(edges[a], x[a]);
            return i;
        }

        static FillTable load(const std::string & path)
        {
            std::ifstream f(path);
            if(!f) throw std::runtime_error("cannot open fill table " + path);
            std::stringstream ss; ss << f.rdbuf();
            const std::string text = ss.str();
            JsonValue j = JsonParser(text).parse();
            FillTable t;
            const JsonValue & e = j.at("edges");
            if(e.kind != JsonValue::kArray || e.arr.size() != 4) throw std::runtime_error("fill table: 'edges' must hold 4 axes");
            for(const auto & ax : e.arr)
            {
                t.edges.push_back(ax.numbers());
                if(t.edges.back().size() < 2) throw std::runtime_error("fill table: an axis has < 2 edges");
            }
            t.frac_contained = j.at("frac_contained").numbers();
            t.frac_exiting = j.at("frac_exiting").numbers();
            if(t.frac_contained.size() != t.ncell() || t.frac_exiting.size() != t.ncell())
                throw std::runtime_error("fill table: fraction size does not match the edges");
            std::vector<double> w = j.at("win_exiting").numbers();
            if(w.size() != 2) throw std::runtime_error("fill table: win_exiting must be [lo, hi]");
            t.win_x_lo = w[0]; t.win_x_hi = w[1];
            if(j.has("win_contained"))
            {
                std::vector<double> c = j.at("win_contained").numbers();
                if(c.size() != 2) throw std::runtime_error("fill table: win_contained must be [lo, hi]");
                t.has_win_c = true; t.win_c_lo = c[0]; t.win_c_hi = c[1];
            }
            if(j.has("elep_min")) t.elep_min = j.at("elep_min").num;
            return t;
        }
    };

    /** Per-event inputs (energies in GeV). */
    struct FillInput
    {
        double containment = NAN;   // leading-muon containment flag (1 = contained)
        double elep = NAN;          // reco E_lep (contained)
        double evis = NAN;          // reco E_vis (contained window, if any)
        double cos = NAN, ehad = NAN, rest = NAN;
        double emu_exiting = NAN;   // exiting E_mu (proxy or MCS), already scaled to GeV
        double vx = NAN, vy = NAN, vz = NAN;
    };

    struct FillResult
    {
        int topology = -1;          // -1 no/invalid muon, 0 contained, 1 exiting
        int window_pass = 0;        // inputs valid and inside the window (contained: always, unless win_contained)
        int pass = 0;               // window_pass && u < f            -> cut option
        double u = -1.0;
        double frac = 0.0;          // cell keep fraction
        double weight = 0.0;        // window_pass ? f : 0              -> weight option
        double elep = NAN, evis = NAN;   // E_lep / E_vis used for the cell and window
    };

    inline FillResult evaluate(const FillTable & t, const FillInput & in)
    {
        FillResult r;
        if(!finite(in.containment)) return r;
        r.topology = (in.containment == 1.0) ? 0 : 1;
        if(finite(in.vx) && finite(in.vy) && finite(in.vz))
            r.u = vertex_uniform(float(in.vx), float(in.vy), float(in.vz));
        const bool contained = (r.topology == 0);
        r.elep = contained ? in.elep : in.emu_exiting;
        r.evis = contained ? in.evis : in.ehad + in.emu_exiting;
        if(!finite(r.elep) || !finite(in.cos) || !finite(in.ehad) || !finite(in.rest) || r.u < 0)
            return r;
        bool win = true;
        if(contained)
        {
            if(t.has_win_c) win = finite(r.evis) && r.evis >= t.win_c_lo && r.evis < t.win_c_hi;
        }
        else
            win = r.elep >= t.elep_min && r.evis >= t.win_x_lo && r.evis < t.win_x_hi;
        const size_t c = t.cell(r.elep, in.cos, in.ehad, in.rest);
        r.frac = contained ? t.frac_contained[c] : t.frac_exiting[c];
        r.window_pass = win ? 1 : 0;
        r.weight = win ? std::max(r.frac, 0.0) : 0.0;   // not clipped at 1: f > 1 is a real weight (the cut can't do it)
        r.pass = (win && r.u < r.frac) ? 1 : 0;
        return r;
    }
}
#endif // GUNDAM_FILL_H
