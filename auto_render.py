# Auto-render PDB structures in a reference style.
# Usage:
#   Save PNGs (headless):   pymol -cq auto_render.py -- --pdbs 5YUI
#   Open in GUI to edit:    pymol auto_render.py -- --pdbs 5YUI --interactive
import argparse
import sys
from pymol import cmd

DATA_DIR = 'data/'
RESULT_DIR = 'results/'

# Reference-style palette
cmd.set_color('helix_purple', [0.60, 0.55, 0.78])
cmd.set_color('sheet_gold',   [0.75, 0.55, 0.30])
cmd.set_color('line_grey',    [0.30, 0.30, 0.30])

# Cel-shaded look: black outlines + rounded helices (visible only in raytraced PNG)
cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_color', 'black')
cmd.set('ray_trace_gain', 0.3)
cmd.set('cartoon_fancy_helices', 1)

parser = argparse.ArgumentParser(description='Render PDB structures in reference style.')
parser.add_argument('--pdbs', nargs='+', default=['1QJG'],
                    help='PDB IDs (without .pdb) to render from DATA_DIR. Example: --pdbs 1QJG 2ABC')
parser.add_argument('--interactive', action='store_true',
                    help='Apply styling and leave the structure loaded in PyMOL (skip ray/PNG/delete).')
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
    cmd.set('line_width', 2.0)

    cmd.orient(pdb)
    cmd.zoom(pdb, buffer=2)

    if args.interactive:
        continue

    cmd.ray(2000, 2000)
    cmd.png(RESULT_DIR + f'{pdb}.png', dpi=300)
    cmd.delete(pdb)
