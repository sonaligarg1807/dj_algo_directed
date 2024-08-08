#another version which reads data from parent bash script directly
import subprocess
import os
import shutil
import time
import sys

# Check for the correct number of command-line arguments
if len(sys.argv) != 3:
    print("Usage: python3 weights_source_resid_dj.py <source_site> <cutoff_distance>")
    sys.exit(1)

# Set variables from command-line arguments
source_site = sys.argv[1]
cutoff_distance = sys.argv[2]

# Set variables for file names and directories
MDP_FILE = "/data/sgarg/test_crystal_plane/5_namd_ham/vector_dj/input_files/namd-qmmm.mdp"
GRO_FILE = "/data/sgarg/test_crystal_plane/5_namd_ham/vector_dj/input_files/nvt.gro"
TOP_FILE = "/data/sgarg/test_crystal_plane/5_namd_ham/vector_dj/input_files/pen-esp.top"
CHARGE_TRANSFER_SCRIPT = "test.sh"
MOLECULE_SPEC = "pentacene-hole.spec"
HAMILTONIAN_FILE = "TB_HAMILTONIAN.xvg"
EXTRACTED_FILE = "extracted.gro"
OUTPUT_FILE = "average_coupling_values.txt"
HAM_LOG_FILE = "ham.log"

# Function to generate the charge_transfer.dat file
def generate_charge_transfer_dat(source_site, target_site):
    charge_transfer_dat = f"""seed = 1
hamiltonian = dftb
tfermi = 300
slkopath = /home/wxie/test-parameters/MIO/
wavefunctionreal = 1.0 0.0
chargecarrier = hole
offdiagscaling = yes
extchrmode = vacuo
espscaling = 1
nsitetypes = 1
sic = 0.0
typefiles = {MOLECULE_SPEC}
nsites = 2
zonesize = 2
sites = {source_site} {target_site}
sitetypes = 1 1
foshift = 0.0 0.0
sitescc = 0 0
jobtype = NOM
deltaqmode = mulliken
internalrelax = no"""
    return charge_transfer_dat

# Function to calculate average coupling from TB_HAMILTONIAN.xvg
def calculate_average_coupling(directory):
    hamiltonian_file = os.path.join(directory, HAMILTONIAN_FILE)
    if os.path.exists(hamiltonian_file):
        with open(hamiltonian_file, 'r') as f:
            lines = f.readlines()
            coupling_values = [float(line.split()[2]) for line in lines if line.strip()]
            avg_coupling = sum(coupling_values) / len(coupling_values) if coupling_values else 0.0
            return abs(avg_coupling)
    else:
        print(f"TB_HAMILTONIAN.xvg file not found in directory: {directory}")
        return None

# Function to check if a job is finished by looking for 'Finished' in ham.log
def is_job_finished(directory):
    log_file = os.path.join(directory, HAM_LOG_FILE)
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            return 'Finished' in f.read()
    return False

# Extract molecules around the source site
subprocess.run(["python3", "extract_molecules_dj.py", GRO_FILE, EXTRACTED_FILE, source_site, cutoff_distance])

# Read the output residues from the extracted file
with open(EXTRACTED_FILE, 'r') as f:
    lines = f.readlines()
    output_resids = set(line.split()[0] for line in lines[2:-1] if line.strip())  # Using a set to store unique residues

# Loop through each unique residue and calculate coupling
for target_site in output_resids:
    print(f"Calculating coupling between sites {source_site} and {target_site}...")

    # Create a directory for the pair
    pair_directory = os.path.join(os.getcwd(), f"{source_site}_{target_site}")
    os.makedirs(pair_directory, exist_ok=True)
    
    # Generate charge_transfer.dat for each pair
    charge_transfer_dat = generate_charge_transfer_dat(source_site, target_site)
    charge_transfer_dat_file = os.path.join(pair_directory, "charge-transfer.dat")
    with open(charge_transfer_dat_file, 'w') as f:
        f.write(charge_transfer_dat)

    # Copy necessary files to the pair directory
    shutil.copy(CHARGE_TRANSFER_SCRIPT, pair_directory)
    shutil.copy(MOLECULE_SPEC, pair_directory)
    shutil.copy(MDP_FILE, pair_directory)

    # Run the GROMACS command to generate ham.tpr in the current directory
    print("Running GROMACS grompp...")
    subprocess.run(["/data/fghalami/gromacs-sh-old_Eik/test_plumed/gromacs-sh-old/COUPLED-DYNAMICS/build-tomas-jan2023/src/kernel/grompp", "-f", MDP_FILE, "-c", GRO_FILE, "-p", TOP_FILE, "-o", 'ham.tpr', "-maxwarn", '1'], cwd=pair_directory)
    # Check if the job has already finished
    if is_job_finished(pair_directory):
        print(f"Job already completed for pair {source_site} and {target_site}. Skipping...")
    else:
        # Run the charge transfer calculation script
        print("Running charge transfer calculation...")
        subprocess.run(['/usr/local/run/ge2011.11/bin/linux-x64/qsub', CHARGE_TRANSFER_SCRIPT], cwd=pair_directory)
        print(f"Calculation initiated for pair {source_site} and {target_site} in directory: {pair_directory}.")

    # Wait for the job to finish
    while not is_job_finished(pair_directory):
        print(f"Job still running for pair {source_site} and {target_site}. Checking again in 10 seconds...")
        time.sleep(10)  # Check again in 0.5 minute

    print(f"Job done for pair {source_site} and {target_site}.")
    
    # Check if the output file already exists, and if not, write the header
    if not os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE, 'w') as f:
            f.write("Source_resid\tTarget_resid\tavg_cpl\n")

    # Calculate average coupling for the pair and save to the output file
    avg_coupling = calculate_average_coupling(pair_directory)
    with open(OUTPUT_FILE, 'a') as f:
        if avg_coupling is not None:
            f.write(f"{source_site}          \t{target_site}          \t{avg_coupling * 1000}\n")
        else:
            f.write(f"{source_site}          \t{target_site}          \tNot available\n")
    #         f.write(f"Average coupling between sites {source_site} and {target_site}: {avg_coupling}\n")
    # else:
    #     with open(OUTPUT_FILE, 'a') as f:
    #         f.write(f"Average coupling between sites {source_site} and {target_site}: Not available (TB_HAMILTONIAN.xvg file missing)\n")

    # Remove test.sh to ensure it's not run again in subsequent iterations
    os.remove(os.path.join(pair_directory, CHARGE_TRANSFER_SCRIPT))

    print(f"Average coupling value saved to {OUTPUT_FILE}.")

print("All calculations completed.")



#........................................................................................................................................
# import subprocess
# import os
# import shutil
# import time

# # Set variables for file names and directories
# MDP_FILE = "/data/sgarg/test_crystal_plane/5_namd_ham/namd/namd-qmmm.mdp"
# GRO_FILE = "/data/sgarg/test_crystal_plane/5_namd_ham/namd/em.gro"
# TOP_FILE = "/data/sgarg/test_crystal_plane/5_namd_ham/namd/pen-esp.top"
# CHARGE_TRANSFER_SCRIPT = "test.sh"
# MOLECULE_SPEC = "pentacene-hole.spec"
# HAMILTONIAN_FILE = "TB_HAMILTONIAN.xvg"
# EXTRACTED_FILE = "extracted.gro"
# OUTPUT_FILE = "average_coupling_values1.txt"
# HAM_LOG_FILE = "ham.log"

# # Function to generate the charge_transfer.dat file
# def generate_charge_transfer_dat(source_site, target_site):
#     charge_transfer_dat = f"""seed = 1
# hamiltonian = dftb
# tfermi = 300
# slkopath = /home/wxie/test-parameters/MIO/
# wavefunctionreal = 1.0 0.0
# chargecarrier = hole
# offdiagscaling = yes
# extchrmode = vacuo
# espscaling = 1
# nsitetypes = 1
# sic = 0.0
# typefiles = {MOLECULE_SPEC}
# nsites = 2
# zonesize = 2
# sites = {source_site} {target_site}
# sitetypes = 1 1
# foshift = 0.0 0.0
# sitescc = 0 0
# jobtype = NOM
# deltaqmode = mulliken
# internalrelax = no"""
#     return charge_transfer_dat

# # Function to calculate average coupling from TB_HAMILTONIAN.xvg
# def calculate_average_coupling(directory):
#     hamiltonian_file = os.path.join(directory, HAMILTONIAN_FILE)
#     if os.path.exists(hamiltonian_file):
#         with open(hamiltonian_file, 'r') as f:
#             lines = f.readlines()
#             coupling_values = [float(line.split()[2]) for line in lines if line.strip()]
#             avg_coupling = sum(coupling_values) / len(coupling_values) if coupling_values else 0.0
#             return avg_coupling
#     else:
#         print(f"TB_HAMILTONIAN.xvg file not found in directory: {directory}")
#         return None

# # Function to check if a job is finished by looking for 'Finished' in ham.log
# def is_job_finished(directory):
#     log_file = os.path.join(directory, HAM_LOG_FILE)
#     if os.path.exists(log_file):
#         with open(log_file, 'r') as f:
#             return 'Finished' in f.read()
#     return False

# # Prompt user for the source site index
# source_site = input("Enter the source site index: ")

# # Extract molecules around the source site
# subprocess.run(["python3", "/data/sgarg/test_crystal_plane/5_namd_ham/namd/extract_molecules.py", GRO_FILE, EXTRACTED_FILE, source_site, "0.25"])

# # Read the output residues from the extracted file
# with open(EXTRACTED_FILE, 'r') as f:
#     lines = f.readlines()
#     output_resids = set(line.split()[0] for line in lines[2:-1] if line.strip())  # Using a set to store unique residues

# # Loop through each unique residue and calculate coupling
# for target_site in output_resids:
#     print(f"Calculating coupling between sites {source_site} and {target_site}...")

#     # Create a directory for the pair
#     pair_directory = os.path.join(os.getcwd(), f"{source_site}_{target_site}")
#     os.makedirs(pair_directory, exist_ok=True)
    
#     # Generate charge_transfer.dat for each pair
#     charge_transfer_dat = generate_charge_transfer_dat(source_site, target_site)
#     charge_transfer_dat_file = os.path.join(pair_directory, "charge-transfer.dat")
#     with open(charge_transfer_dat_file, 'w') as f:
#         f.write(charge_transfer_dat)

#     # Copy necessary files to the pair directory
#     shutil.copy(CHARGE_TRANSFER_SCRIPT, pair_directory)
#     shutil.copy(MOLECULE_SPEC, pair_directory)
#     shutil.copy(MDP_FILE, pair_directory)

#     # Run the GROMACS command to generate ham.tpr in the current directory
#     print("Running GROMACS grompp...")
#     subprocess.run(["/data/fghalami/gromacs-sh-old_Eik/test_plumed/gromacs-sh-old/COUPLED-DYNAMICS/build-tomas-jan2023/src/kernel/grompp", "-f", MDP_FILE, "-c", GRO_FILE, "-p", TOP_FILE, "-o", 'ham.tpr', "-maxwarn", '1'], cwd=pair_directory)

#     # Check if the job has already finished
#     if is_job_finished(pair_directory):
#         print(f"Job already completed for pair {source_site} and {target_site}. Skipping...")
#     else:
#         # Run the charge transfer calculation script
#         print("Running charge transfer calculation...")
#         subprocess.run(["qsub", CHARGE_TRANSFER_SCRIPT], cwd=pair_directory)
#         print(f"Calculation initiated for pair {source_site} and {target_site} in directory: {pair_directory}.")

#     # Wait for the job to finish
#     while not is_job_finished(pair_directory):
#         print(f"Job still running for pair {source_site} and {target_site}. Checking again in 1 minute...")
#         time.sleep(60)  # Check again in 1 minute

#     print(f"Job done for pair {source_site} and {target_site}.")
    
#     # Check if the output file already exists, and if not, write the header
# if not os.path.exists(OUTPUT_FILE):
#     with open(OUTPUT_FILE, 'w') as f:
#         f.write("Source_resid\tTarget_resid\tavg_cpl\n")

#     # Calculate average coupling for the pair and save to the output file
#     avg_coupling = calculate_average_coupling(pair_directory)
#     with open(OUTPUT_FILE, 'a') as f:
#         if avg_coupling is not None:
#             f.write(f"{source_site}          \t{target_site}          \t{avg_coupling * 1000}\n")
#         else:
#             f.write(f"{source_site}          \t{target_site}          \tNot available\n")
#     #         f.write(f"Average coupling between sites {source_site} and {target_site}: {avg_coupling}\n")
#     # else:
#     #     with open(OUTPUT_FILE, 'a') as f:
#     #         f.write(f"Average coupling between sites {source_site} and {target_site}: Not available (TB_HAMILTONIAN.xvg file missing)\n")

#     # Remove test.sh to ensure it's not run again in subsequent iterations
#     os.remove(os.path.join(pair_directory, CHARGE_TRANSFER_SCRIPT))

#     print(f"Average coupling value saved to {OUTPUT_FILE}.")

# print("All calculations completed.")
