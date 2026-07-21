# setup.sh

# sl7
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh

# For ICARUS:
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

# Set up the required dependencies:
setup sbnana v10_01_02_01 -q e26:prof

setup cmake v3_27_4

# ICARUS
htgettoken -a htvaultprod.fnal.gov  -i icarus

# Running in the grid

# Create Project
rm -fr /pnfs/icarus/scratch/users/faabdalr/New_med/numu_t/
python3 batch/medulla.py -t /exp/icarus/data/users/fatima/medulla_numu/selection/toml/numuCC_inclusive_uncontained.toml -p /pnfs/icarus/scratch/users/faabdalr/New_med/numu_t -b 1 --create-project

# Submit a test
python3 batch/medulla.py -p /pnfs/icarus/scratch/users/faabdalr/New_med/numu_t -e icarus --test-job --branch feature/numu_analysis --memory 12000 --disk 30 --lifetime 2

# Submit all jobs
python3 batch/medulla.py -p /pnfs/icarus/scratch/users/faabdalr/New_med/numu_f -e icarus --launch-jobs --branch feature/numu_analysis --memory 12000 --disk 30 --lifetime 10

# Verify the test output 
ls -lh /pnfs/icarus/scratch/users/faabdalr/nueCC_inclusive_all_ntest/output/

#Open the ROOT file and check for TTrees with entries:
# spack env
. /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh

# root
# - hljl7gy root@6.28.12%gcc@12.2.0 arch=linux-almalinux9-x86_64_v3
spack load /hljl7gy

root -l /nueCC_inclusive_ext_ntest/output/output_jobid0000.root

#Inside ROOT:
_file0->ls()
events->cd("NuMIFull")
selected->GetEntries()


# Monitor 
jobsub_q -G icarus --jobid <JOBID>
jobsub_q -G icarus --jobid 85359048.0@jobsub01.fnal.gov

condor_q 70076075 -name jobsub03.fnal.gov


# Merge output files
hadd -f /exp/icarus/data/users/faabdalr/medulla_output/NueCC_Inclusive.root /pnfs/icarus/scratch/users/faabdalr/NueCC_Inclusive/output/*.root



# Clone the medulla repository:
git clone https://github.com/justinjmueller/medulla.git medulla
cd medulla && git checkout v1.0.3
mkdir build && cd build

# Configure and build medulla:
setup cmake v3_27_4
cmake .. && make -j4


# Running the Selection
./selection/medulla /exp/icarus/data/users/fatima/medulla_numu/selection/toml/nu

