# Interactive view setter: load PDBs with the auto_render styling, rotate
# manually in the GUI, then save the orientation for later replay.
#
# Usage:
#   pymol set_view.py -- --pdbs 5AOU [--motif]
#
# In the PyMOL command line:
#   view_pdb <PDB>    Switch which loaded PDB is visible / oriented (when multiple).
#   save_view <PDB>   Save the current view to views/<PDB>.txt. (e.g. save_view 5AOU)
import argparse
import os
import sys
from pymol import cmd

from utils import default_task_name, get_motif_info, is_benchmark_pdb

DATA_DIR = 'data/'
VIEWS_DIR = 'views/'

cmd.set_color('helix_purple', [0.60, 0.55, 0.78])
cmd.set_color('sheet_gold',   [0.75, 0.55, 0.30])
cmd.set_color('line_grey',    [0.33, 0.33, 0.33])
cmd.set_color('motif_red',    [0.80, 0.40, 0.40])

cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_color', 'black')
cmd.set('ray_trace_gain', 0.3)
cmd.set('cartoon_fancy_helices', 1)

parser = argparse.ArgumentParser(description='Set and save custom views for PDB structures.')
parser.add_argument('--pdbs', nargs='+', required=True,
                    help='PDB IDs (without .pdb) to load from DATA_DIR.')
parser.add_argument('--motif', action='store_true',
                    help='Highlight motif residues for registered benchmark PDBs.')
args = parser.parse_args(sys.argv[1:])
pdb_list = args.pdbs

os.makedirs(VIEWS_DIR, exist_ok=True)


def load_saved_view(pdb):
    """Return the 18-float view tuple from views/<pdb>.txt, or None if missing."""
    path = os.path.join(VIEWS_DIR, f'{pdb}.txt')
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return tuple(float(line) for line in f if line.strip())


def apply_default_view(pdb):
    """Use the saved view if available, otherwise the default orient + zoom."""
    saved = load_saved_view(pdb)
    if saved is not None:
        cmd.set_view(saved)
        print(f'[set_view] Loaded saved view for {pdb} from {VIEWS_DIR}{pdb}.txt')
    else:
        cmd.orient(pdb)
        cmd.zoom(pdb, buffer=2)


for pdb in pdb_list:
    cmd.load(DATA_DIR + pdb + '.pdb', pdb)

    cmd.bg_color('white')
    cmd.hide('everything', pdb)
    cmd.dss(pdb)

    cmd.show('cartoon', pdb)
    cmd.set('cartoon_color', 'helix_purple', pdb)
    cmd.set('cartoon_color', 'sheet_gold', f'{pdb} and ss S')

    cmd.show('lines', pdb)
    cmd.color('line_grey', pdb)
    cmd.set('line_width', 1.6)

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

    if pdb != pdb_list[0]:
        cmd.disable(pdb)

apply_default_view(pdb_list[0])


def view_pdb(pdb):
    """Show only `pdb`, loading its saved view if present, else default orient + zoom."""
    if pdb not in pdb_list:
        print(f"[set_view] '{pdb}' is not loaded. Loaded: {', '.join(pdb_list)}")
        return
    for p in pdb_list:
        cmd.disable(p)
    cmd.enable(pdb)
    apply_default_view(pdb)


def save_view(pdb):
    """Save the current view to views/<pdb>.txt as one float per line."""
    if pdb not in pdb_list:
        print(f"[set_view] '{pdb}' is not loaded. Loaded: {', '.join(pdb_list)}")
        return
    view = cmd.get_view()
    path = os.path.join(VIEWS_DIR, f'{pdb}.txt')
    with open(path, 'w') as f:
        for v in view:
            f.write(f'{v}\n')
    print(f'[set_view] Saved current view to {path}')


cmd.extend('view_pdb', view_pdb)
cmd.extend('save_view', save_view)

print('[set_view] Rotate the structure with the mouse, then in the PyMOL command line:')
print('[set_view]   save_view <PDB>      e.g. save_view ' + pdb_list[0])
if len(pdb_list) > 1:
    print('[set_view]   view_pdb <PDB>       switch to another loaded PDB')
print(f'[set_view] Loaded PDBs: {", ".join(pdb_list)}')
