# MM-GBSA workflow (step-by-step) for HSA + 24 ligands
# Put this code in a Jupyter notebook cell-by-cell. No bash cells used.
# IMPORTANT: This notebook calls external AmberTools executables (antechamber, parmchk2, tleap, MMPBSA.py, sander/pmemd).
# You must have AmberTools installed (conda install -c conda-forge -c ambermd ambertools) for those commands to work.

# %%
"""
Cell 1: Imports and user paths (paste your filepaths here)
"""
import os
import sys
import shutil
import subprocess
from pathlib import Path
import textwrap
import json

# --- Proteins (as provided by user) ---
h2bxd = "HSA_Target/2BXD/2bxd.pdb"
h2bxg = "HSA_Target/2BXD/2bxg.pdb"
h1ao6 = "HSA_Target/1AO6/1ao6.pdb"

p2bxd = "HSA_Target/2BXD/prep_2bxd.pdb"
p2bxg = "HSA_Target/2BXD/prep_2bxg.pdb"
p1ao6 = "HSA_Target/1AO6/prep_1ao6.pdb"

# --- Ligands (as provided by user) ---
# p2bxd (Site1)
h1ea = "Results/Site1/HOLO/ENA/ENA-Site1_1h.pdb"
h1ea1 = "Results/Site1/HOLO/ENA1/ENA1-Site1_1h.pdb"
h1eb = "Results/Site1/HOLO/ENB/ENB-Site1_1h.pdb"
h1eb1 = "Results/Site1/HOLO/ENB1/ENB1-Site1_1h.pdb"
h1rwf = "Results/Site1/HOLO/RWF/RWF-Site1_1h.pdb"
h1rbp = "Results/Site1/HOLO/IBP/IBP-Site1_1h.pdb"

# p2bxg (Site2)
h2ea = "Results/Site2/HOLO/ENA/ENA-Site2_1h.pdb"
h2ea1 = "Results/Site2/HOLO/ENA1/ENA1-Site2_1h.pdb"
h2eb = "Results/Site2/HOLO/ENB/ENB-Site2_1h.pdb"
h2eb1 = "Results/Site2/HOLO/ENB1/ENB1-Site2_1h.pdb"
h2rwf = "Results/Site2/HOLO/RWF/RWF-Site2_1h.pdb"
h2ibp = "Results/Site2/HOLO/IBP/IBP-Site2_1h.pdb"

# p1ao6 (APO Site1)
a1ea = "Results/Site1/APO/ENA/ENA-Site1_1h.pdb"
a1ea1 = "Results/Site1/APO/ENA1/ENA1-Site1_1h.pdb"
a1eb = "Results/Site1/APO/ENB/ENB-Site1_1h.pdb"
a1eb1 = "Results/Site1/APO/ENB1/ENB1-Site1_1h.pdb"
a1rwf = "Results/Site1/APO/RWF/RWF-Site1_1h.pdb"
a1ibp = "Results/Site1/APO/IBP/IBP-Site1_1h.pdb"

# p1ao6 (APO Site2)
a2ea = "Results/Site2/APO/ENA/ENA-Site2_1h.pdb"
a2ea1 = "Results/Site2/APO/ENA1/ENA1-Site2_1h.pdb"
a2eb = "Results/Site2/APO/ENB/ENB-Site2_1h.pdb"
a2eb1 = "Results/Site2/APO/ENB1/ENB1-Site2_1h.pdb"
a2rwf = "Results/Site2/APO/RWF/RWF-Site2_1h.pdb"
a2ibp = "Results/Site2/APO/IBP/IBP-Site2_1h.pdb"

# results dir
results_dir = Path("Results/MM")
results_dir.mkdir(parents=True, exist_ok=True)

# Build a list of (protein_prep, ligand_path, label) tuples to iterate
jobs = [
    # using p2bxd
    (p2bxd, h1ea, "p2bxd_h1ea"),
    (p2bxd, h1ea1, "p2bxd_h1ea1"),
    (p2bxd, h1eb, "p2bxd_h1eb"),
    (p2bxd, h1eb1, "p2bxd_h1eb1"),
    (p2bxd, h1rwf, "p2bxd_h1rwf"),
    (p2bxd, h1rbp, "p2bxd_h1rbp"),
    # p2bxg
    (p2bxg, h2ea, "p2bxg_h2ea"),
    (p2bxg, h2ea1, "p2bxg_h2ea1"),
    (p2bxg, h2eb, "p2bxg_h2eb"),
    (p2bxg, h2eb1, "p2bxg_h2eb1"),
    (p2bxg, h2rwf, "p2bxg_h2rwf"),
    (p2bxg, h2ibp, "p2bxg_h2ibp"),
    # p1ao6 site1
    (p1ao6, a1ea, "p1ao6_a1ea"),
    (p1ao6, a1ea1, "p1ao6_a1ea1"),
    (p1ao6, a1eb, "p1ao6_a1eb"),
    (p1ao6, a1eb1, "p1ao6_a1eb1"),
    (p1ao6, a1rwf, "p1ao6_a1rwf"),
    (p1ao6, a1ibp, "p1ao6_a1ibp"),
    # p1ao6 site2
    (p1ao6, a2ea, "p1ao6_a2ea"),
    (p1ao6, a2ea1, "p1ao6_a2ea1"),
    (p1ao6, a2eb, "p1ao6_a2eb"),
    (p1ao6, a2eb1, "p1ao6_a2eb1"),
    (p1ao6, a2rwf, "p1ao6_a2rwf"),
    (p1ao6, a2ibp, "p1ao6_a2ibp"),
]

print(f"Prepared {len(jobs)} jobs. Results will be written to: {results_dir.resolve()}")

# %%
"""
Cell 2: Helper utilities - check for executables, run commands, and (optionally) pip install Python-only dependencies.
"""
import shutil

def which(program):
    """Like UNIX which"""
    return shutil.which(program)


def run_cmd(cmd, check=True):
    print("RUN:", " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.stdout:
        print(proc.stdout)
    if proc.stderr:
        print(proc.stderr)
    if check and proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nReturn code: {proc.returncode}")
    return proc


# Check for AmberTools executables (antechamber, parmchk2, tleap, MMPBSA.py)
required_bins = ["antechamber", "parmchk2", "tleap", "MMPBSA.py"]
found = {b: which(b) for b in required_bins}
print(json.dumps(found, indent=2))

# Attempt to pip install python-only dependencies (pytraj, parmed, mdtraj, biopython)
# These are optional helpers for analysis. AmberTools still required separately.
python_pkgs = ["parmed", "mdtraj", "pytraj", "biopython", "numpy", "pandas"]
print('\nYou can install Python helper packages (dry-run below).')
print('If you want to auto-install, uncomment the pip-install block in the cell.')

# If desired, uncomment below to auto-install via pip from the notebook (not recommended for ambertools):
# for pkg in python_pkgs:
#     run_cmd([sys.executable, "-m", "pip", "install", pkg])

# AmberTools note
if not found['tleap']:
    print('\n\nNOTE: tleap not found in PATH. You need AmberTools (conda install -c conda-forge -c ambermd ambertools).')
    print('If you have conda, run (in terminal):\n  conda install -c conda-forge -c ambermd ambertools')

# Check for pmemd.cuda or sander (GPU vs CPU MD)
pmemd_cuda = which('pmemd.cuda')
sander_bin = which('sander')
print(f"pmemd.cuda: {pmemd_cuda}, sander: {sander_bin}")

# %%
"""
Cell 3: Utility functions to parameterize ligand (antechamber), build topologies with tleap, and optionally run quick single-frame MM-GBSA.
"""

def parametrize_ligand(ligand_pdb, out_dir):
    """Runs antechamber + parmchk2 to generate mol2 + frcmod. Returns paths.
    ligand_pdb: path to ligand coordinates (PDB) - preferably original ligand with proper connectivity.
    out_dir: directory to write files.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = Path(ligand_pdb).stem
    mol2 = out_dir / f"{base}.mol2"
    frcmod = out_dir / f"{base}.frcmod"

    # Step 1: antechamber
    cmd1 = ["antechamber", "-i", str(ligand_pdb), "-fi", "pdb", "-o", str(mol2), "-fo", "mol2", "-c", "bcc", "-s", "2"]
    run_cmd(cmd1)

    # Step 2: parmchk2
    cmd2 = ["parmchk2", "-i", str(mol2), "-f", "mol2", "-o", str(frcmod)]
    run_cmd(cmd2)

    return str(mol2), str(frcmod)


def write_tleap_and_run(protein_pdb, ligand_mol2, ligand_frcmod, out_prefix, solvate=True, box=10.0):
    """Write tleap input that creates receptor, ligand, complex prmtops and solvated complex.
    out_prefix is the prefix for saving files.
    Returns paths to complex.prmtop, receptor.prmtop, ligand.prmtop, complex_solvated.prmtop
    """
    tleap_in = f"""
source leaprc.protein.ff14SB
source leaprc.gaff

rec = loadpdb {protein_pdb}
lig = loadmol2 {ligand_mol2}
loadamberparams {ligand_frcmod}

complex = combine {{rec lig}}

saveamberparm rec {out_prefix}_receptor.prmtop {out_prefix}_receptor.inpcrd
saveamberparm lig {out_prefix}_ligand.prmtop {out_prefix}_ligand.inpcrd
saveamberparm complex {out_prefix}_complex.prmtop {out_prefix}_complex.inpcrd

"""
    if solvate:
        tleap_in += f"solvatebox complex TIP3PBOX {box}\naddions complex Na+ 0\nsaveamberparm complex {out_prefix}_complex_solvated.prmtop {out_prefix}_complex_solvated.inpcrd\n"
    tleap_in += "quit\n"

    tleap_file = Path(out_prefix + "_tleap.in")
    tleap_file.write_text(tleap_in)
    run_cmd(["tleap", "-f", str(tleap_file)])

    # return prmtop paths (strings)
    cp = f"{out_prefix}_complex.prmtop"
    rp = f"{out_prefix}_receptor.prmtop"
    lp = f"{out_prefix}_ligand.prmtop"
    sp = f"{out_prefix}_complex_solvated.prmtop" if solvate else None
    return cp, rp, lp, sp


def run_single_frame_mmgbsa(cp, rp, lp, out_prefix, mmpbsa_in=None):
    """Run single-frame MM-GBSA (no -y trajectory). Requires MMPBSA.py in PATH.
    Returns path to output file.
    """
    if mmpbsa_in is None:
        mmpbsa_in = 'mmpbsa.in'
        Path(mmpbsa_in).write_text(textwrap.dedent('''
&general
  verbose=1,
/
&gb
  igb=5, saltcon=0.150
/
        '''))
    out_file = f"{out_prefix}_MMPBSA.dat"
    cmd = ["MMPBSA.py", "-O", "-i", mmpbsa_in, "-o", out_file, "-cp", cp, "-rp", rp, "-lp", lp]
    run_cmd(cmd)
    return out_file

# %%
"""
Cell 4: Example single-job run (fast pathway)
This does:
 - parametrize ligand with antechamber
 - build prmtops with tleap
 - run single-frame MM-GBSA (fast test)

Use this first to screen all 24 ligands quickly before running MD-based MM-GBSA.
"""

def single_job_quick(job, work_root="mm_work"):
    protein_prep, ligand_pdb, label = job
    wd = Path(work_root) / label
    wd.mkdir(parents=True, exist_ok=True)

    # parametrize ligand
    print(f"\n=== JOB: {label} ===")
    mol2, frcmod = parametrize_ligand(ligand_pdb, out_dir=wd)

    # build prmtops
    out_prefix = wd / label
    cp, rp, lp, sp = write_tleap_and_run(protein_prep, mol2, frcmod, str(out_prefix))

    # run single-frame mmgbsa (fast)
    mm_out = run_single_frame_mmgbsa(cp, rp, lp, out_prefix)

    # copy results
    dest = results_dir / f"{label}_singleframe_MMPBSA.dat"
    shutil.copy(mm_out, dest)
    print(f"Saved: {dest}")
    return dest

# Example: run only one job as a test
# test_output = single_job_quick(jobs[0])
# print(test_output)

# %%
"""
Cell 5: Batch-run quick single-frame MM-GBSA for all jobs (uncomment to run)
WARNING: This is CPU-light but still calls AmberTools for each ligand. It is a good screening step.
"""

# Uncomment to run batch quick-mode across all jobs
# for job in jobs:
#     try:
#         single_job_quick(job)
#     except Exception as e:
#         print(f"Job {job[2]} failed: {e}")

# %%
"""
Cell 6: If you want MD-based MM-GBSA (recommended for highest quality) - minimal pipeline outline.
This code *writes* the required input files and shows commands to run. Running full MD for 24 ligands is expensive.

The functions below create MD input files and call pmemd.cuda (if available) or sander (CPU).
"""

def write_min_heat_npt_prod(wd, prefix):
    wd = Path(wd)
    min_in = textwrap.dedent('''
Minimization
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=8.0
 /
''')
    heat_in = textwrap.dedent('''
Heating
 &cntrl
  imin=0, irest=0, ntx=1,
  ntt=3, temp0=300.0, tempi=0.0,
  nstlim=50000, dt=0.002, ntpr=500, ntwx=500,
  ntb=1, cut=8.0
 /
''')
    npt_in = textwrap.dedent('''
NPT equil
 &cntrl
  imin=0, irest=1, ntx=5,
  ntt=3, temp0=300.0,
  nstlim=50000, dt=0.002, ntpr=500, ntwx=500,
  ntb=2, pres0=1.0, taup=2.0, cut=8.0
 /
''')
    prod_in = textwrap.dedent('''
Production
 &cntrl
  imin=0, irest=1, ntx=5,
  ntt=3, temp0=300.0,
  nstlim=500000, dt=0.002, ntpr=1000, ntwx=1000,
  ntb=2, pres0=1.0, taup=2.0, cut=8.0
 /
''')
    (wd / f"{prefix}_min.in").write_text(min_in)
    (wd / f"{prefix}_heat.in").write_text(heat_in)
    (wd / f"{prefix}_npt.in").write_text(npt_in)
    (wd / f"{prefix}_prod.in").write_text(prod_in)
    return [wd / f"{prefix}_min.in", wd / f"{prefix}_heat.in", wd / f"{prefix}_npt.in", wd / f"{prefix}_prod.in"]


def run_md_pipeline(wd, prefix, prmtop, inpcrd):
    """Run MD stages. Detect pmemd.cuda -> use it, else use sander.
    wd: working dir
    prmtop: solvated prmtop
    inpcrd: initial inpcrd
    """
    wd = Path(wd)
    pmemd = which('pmemd.cuda') or which('pmemd')
    sander = which('sander')
    mdbin = pmemd or sander
    if not mdbin:
        raise RuntimeError('No MD engine found (pmemd.cuda or sander). Install Amber or use GROMACS alternative).')

    prefix = prefix
    # run minimize
    run_cmd([mdbin, "-O", "-i", str(wd / f"{prefix}_min.in"), "-o", str(wd / f"{prefix}_min.out"), "-p", prmtop, "-c", inpcrd, "-r", str(wd / f"{prefix}_min.rst"), "-ref", inpcrd])
    # heat
    run_cmd([mdbin, "-O", "-i", str(wd / f"{prefix}_heat.in"), "-o", str(wd / f"{prefix}_heat.out"), "-p", prmtop, "-c", str(wd / f"{prefix}_min.rst"), "-r", str(wd / f"{prefix}_heat.rst"), "-x", str(wd / f"{prefix}_heat.nc"), "-ref", inpcrd])
    # npt
    run_cmd([mdbin, "-O", "-i", str(wd / f"{prefix}_npt.in"), "-o", str(wd / f"{prefix}_npt.out"), "-p", prmtop, "-c", str(wd / f"{prefix}_heat.rst"), "-r", str(wd / f"{prefix}_npt.rst"), "-x", str(wd / f"{prefix}_npt.nc")])
    # prod
    run_cmd([mdbin, "-O", "-i", str(wd / f"{prefix}_prod.in"), "-o", str(wd / f"{prefix}_prod.out"), "-p", prmtop, "-c", str(wd / f"{prefix}_npt.rst"), "-r", str(wd / f"{prefix}_prod.rst"), "-x", str(wd / f"{prefix}_prod.nc")])
    return str(wd / f"{prefix}_prod.nc")

# %%
"""
Cell 7: Example pipeline that does MD-based MM-GBSA for ONE selected job (manual step-by-step run recommended).
Run this only after testing the quick single-frame pipeline and after you confirm Amber MD executables are present.
"""

def md_job_example(job, work_root='mm_work_md'):
    protein_prep, ligand_pdb, label = job
    wd = Path(work_root) / label
    wd.mkdir(parents=True, exist_ok=True)

    # 1) parametrize ligand
    mol2, frcmod = parametrize_ligand(ligand_pdb, wd)

    # 2) tleap to create solvated system
    out_prefix = wd / label
    cp, rp, lp, sp = write_tleap_and_run(protein_prep, mol2, frcmod, str(out_prefix), solvate=True, box=10.0)

    # 3) write md inputs
    write_min_heat_npt_prod(wd, label)

    # 4) run MD pipeline (will use pmemd.cuda if available)
    traj = run_md_pipeline(wd, label, sp, str(wd / f"{label}_complex_solvated.inpcrd"))

    # 5) run MMPBSA on trajectory
    mmpbsa_in = 'mmpbsa_md.in'
    Path(mmpbsa_in).write_text(textwrap.dedent('''
&general
  startframe=1, endframe=50, interval=5,
  verbose=1,
/
&gb
  igb=5, saltcon=0.150
/
    '''))
    out_file = f"{label}_md_MMPBSA.dat"
    run_cmd(["MMPBSA.py", "-O", "-i", mmpbsa_in, "-o", out_file, "-sp", sp, "-cp", cp, "-rp", rp, "-lp", lp, "-y", traj])

    # copy results
    dest = results_dir / out_file
    shutil.copy(out_file, dest)
    print(f"MD MMPBSA saved to: {dest}")
    return dest

# Don't auto-run. Use md_job_example(jobs[0]) to test a single job.

# %%
"""
Cell 8: Postprocessing - collect all single-frame outputs and summarize.
"""

def collect_results(results_dir=results_dir):
    results_dir = Path(results_dir)
    summary = []
    for path in results_dir.glob('*MMPBSA*.dat'):
        text = path.read_text()
        # crude parsing: find lines with "DELTA TOTAL" or "TOTAL" depending on output format
        for line in text.splitlines():
            if 'DELTA TOTAL' in line or line.strip().startswith('TOTAL'):
                summary.append((path.name, line.strip()))
                break
    return summary

# Example: print collected results
# print(collect_results())

# %%
"""
Cell 9: Notes & troubleshooting (short).
- If you get OpenBabel kekulization warnings when converting .pdbqt -> .pdb/.mol2, regenerate ligand from original source (SDF/MOL2) if possible.
- For parameterization, antechamber assumes reasonable chemistry; check atom types in the generated mol2 before tleap.
- For MD speed on GPU: pmemd.cuda (Amber) is fastest but requires Amber package and a license for pmemd.cuda in some configurations. pmemd.cuda might be present if you have Amber with GPU support.
- Single-frame MM-GBSA (quick) is good for screening 24 ligands. MD-based MM-GBSA is more reliable but costly.

# End of notebook script.
