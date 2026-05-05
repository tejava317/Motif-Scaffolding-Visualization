from pymol import cmd

DATA_DIR = 'data/'
RESULT_DIR = 'results/'

pdb_list = ['1QJG']

for pdb in pdb_list:
    # Load the PDB file
    cmd.load(DATA_DIR + pdb + '.pdb')

    # Process the structure
    cmd.bg_color('white')
    cmd.orient()
    
    # Render the result
    cmd.ray(2000, 2000)
    cmd.png(RESULT_DIR + f"{pdb}.png")

    # Clean up for the next iteration
    cmd.delete(pdb)
