# Auto-render PDB structures in a reference style.
#
# Usage:
#   Save PNGs (headless):   pymol -cq auto_render.py -- [OPTIONS]
#   Open in GUI to edit:    pymol    auto_render.py -- [OPTIONS]
#
# Options:
#   --pdbs ID [ID ...]   PDB IDs (without .pdb) to render from data/ (default: 5YUI)
#   --interactive        Skip ray/PNG/delete and leave the styled structure in the GUI
#   --motif              Highlight motif residues in red for registered benchmark PDBs
import argparse
import sys
from pymol import cmd

from utils import default_task_name, get_motif_info, is_benchmark_pdb

DATA_DIR = 'data/'
RESULT_DIR = 'results/'

# Reference-style palette
cmd.set_color('helix_purple', [0.60, 0.55, 0.78])
cmd.set_color('sheet_gold',   [0.75, 0.55, 0.30])
cmd.set_color('line_grey',    [0.33, 0.33, 0.33])
cmd.set_color('motif_red',    [0.80, 0.40, 0.40])

# Cel-shaded look: black outlines + rounded helices (visible only in raytraced PNG)
cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_color', 'black')
cmd.set('ray_trace_gain', 0.3)
cmd.set('cartoon_fancy_helices', 1)

parser = argparse.ArgumentParser(description='Render PDB structures in reference style.')
parser.add_argument('--pdbs', nargs='+', default=['5YUI'],
                    help='PDB IDs (without .pdb) to render from DATA_DIR. Example: --pdbs 5YUI 2ABC')
parser.add_argument('--interactive', action='store_true',
                    help='Apply styling and leave the structure loaded in PyMOL (skip ray/PNG/delete).')
parser.add_argument('--motif', action='store_true',
                    help='Highlight motif residues in red for registered benchmark PDBs.')
args = parser.parse_args(sys.argv[1:])
pdb_list = args.pdbs

for pdb in pdb_list:
    cmd.load(DATA_DIR + pdb + '.pdb', pdb)

    cmd.bg_color('white')
    cmd.hide('everything', pdb)

    # Recompute secondary structure (this PDB has no HELIX/SHEET records)
    cmd.dss(pdb)

    # Cartoon for whole structure: helices purple, sheets gold
    cmd.show('cartoon', pdb)
    cmd.set('cartoon_color', 'helix_purple', pdb)
    cmd.set('cartoon_color', 'sheet_gold', f'{pdb} and ss S')

    # Thin grey atomic lines — the wire-through-structure look from the references
    cmd.show('lines', pdb)
    cmd.color('line_grey', pdb)
    cmd.set('line_width', 1.6)

    # Motif highlighting: thick sticks on motif residue side chains, red color
    # only on the conditioning atoms (tip atoms) — matches the paper's figure style.
    orient_target = pdb
    if args.motif and is_benchmark_pdb(pdb):
        info = get_motif_info(
            generated_pdb_path=DATA_DIR + pdb + '.pdb',
            task_name=default_task_name(pdb),
        )
        resi_sel = f"{pdb} and resi {info['selection']}"
        atom_sel = f"{pdb} and ({info['atom_selection']})"
        cmd.show('sticks', resi_sel)
        cmd.set('stick_radius', 0.30, resi_sel)
        cmd.color('motif_red', atom_sel)
        # orient_target = resi_sel

    cmd.orient(orient_target)
    cmd.zoom(pdb, buffer=2)

    if args.interactive:
        continue

    cmd.ray(2000, 2000)
    cmd.png(RESULT_DIR + f'{pdb}.png', dpi=300)
    cmd.delete(pdb)
