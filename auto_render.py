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
#   --view               Use the saved orientation in views/<PDB>.txt (created by set_view.py)
import argparse
import os
import sys
from pymol import cmd

from utils import default_task_name, get_motif_info, is_benchmark_pdb

DATA_DIR = 'data/'
RESULT_DIR = 'results/'
MOTIFS_DIR = 'results/motifs/'
VIEWS_DIR = 'views/'

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
parser.add_argument('--pdbs', nargs='+', required=True,
                    help='PDB IDs (without .pdb) to render from DATA_DIR. Example: --pdbs 5YUI 2ABC')
parser.add_argument('--motif', action='store_true',
                    help='Highlight motif residues in red for registered benchmark PDBs.')
parser.add_argument('--interactive', action='store_true',
                    help='Apply styling and leave the structure loaded in PyMOL (skip ray/PNG/delete).')
parser.add_argument('--view', action='store_true',
                    help='Use the saved orientation in views/<PDB>.txt instead of the default orient.')
args = parser.parse_args(sys.argv[1:])
pdb_list = args.pdbs


def load_saved_view(pdb):
    """Return the 18-float view tuple from views/<pdb>.txt, or None if missing."""
    path = os.path.join(VIEWS_DIR, f'{pdb}.txt')
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return tuple(float(line) for line in f if line.strip())

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
    motif_info = None
    if args.motif and is_benchmark_pdb(pdb):
        motif_info = get_motif_info(
            generated_pdb_path=DATA_DIR + pdb + '.pdb',
            task_name=default_task_name(pdb),
        )
        resi_sel = f"{pdb} and resi {motif_info['selection']}"
        atom_sel = f"{pdb} and ({motif_info['atom_selection']})"
        cmd.show('sticks', resi_sel)
        cmd.set('stick_radius', 0.30, resi_sel)
        cmd.color('motif_red', atom_sel)
        # orient_target = resi_sel

    saved_view = load_saved_view(pdb) if args.view else None
    if saved_view is not None:
        cmd.set_view(saved_view)
        print(f'[auto_render] Loaded saved view for {pdb} from {VIEWS_DIR}{pdb}.txt')
    else:
        if args.view:
            print(f'[auto_render] No saved view for {pdb}; falling back to default orient.')
        cmd.orient(orient_target)
        cmd.zoom(pdb, buffer=2)

    if args.interactive:
        print(f'[auto_render] Loaded {pdb} in interactive mode (skipping ray/PNG).')
        continue

    print(f'[auto_render] Rendering {pdb}...')
    cmd.ray(2000, 2000)
    out_path = RESULT_DIR + f'{pdb}.png'
    cmd.png(out_path, dpi=300)
    print(f'[auto_render] Saved {out_path}')

    if motif_info is not None:
        os.makedirs(MOTIFS_DIR, exist_ok=True)
        # Kill depth-cue fog and use ortho projection so the close-up background
        # stays crisp instead of fading toward the (white) background color.
        cmd.set('depth_cue', 0)
        cmd.set('fog', 0)
        cmd.set('ray_trace_fog', 0)
        cmd.set('orthoscopic', 1)
        for idx, resi in enumerate(motif_info['motif_resi']):
            letter = chr(ord('a') + idx)
            cmd.zoom(f"{pdb} and resi {resi}", buffer=3)
            # cmd.zoom tightens the clip slab around the residue, which crops
            # the rest of the protein. Widen it so back helices stay visible.
            cmd.clip('slab', 200)
            print(f'[auto_render] Rendering motif {letter} (resi {resi}) for {pdb}...')
            cmd.ray(2000, 2000)
            motif_out = os.path.join(MOTIFS_DIR, f'{pdb}_motif_{letter}.png')
            cmd.png(motif_out, dpi=300)
            print(f'[auto_render] Saved {motif_out}')
        cmd.set('depth_cue', 1)
        cmd.set('fog', 1)
        cmd.set('ray_trace_fog', 1)
        cmd.set('orthoscopic', 0)

    cmd.delete(pdb)
