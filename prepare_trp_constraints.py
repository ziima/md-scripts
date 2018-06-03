"""Generate a constraints for TRP channel simulations."""
from pyvmd.atoms import Selection
from pyvmd.molecules import Molecule, FORMAT_PDB

PSF = 'initial.psf'
PDB = '00_minimize/00_minimize.coor'
# Lipid tail simulation constraint
LIPID_TAILS_PDB = 'tails_only.pdb'
# Protein fixed simulation constraint
PROTEIN_PDB = 'protein.pdb'
# Protein CA constraint
BACKBONE_CA_PDB = 'backbone_ca.pdb'

model = Molecule.create()
model.load(PSF)
model.load(PDB, FORMAT_PDB)
model_sel = Selection('all', model)
model_sel.atomsel.set('beta', 0)

# Lipid constraints
# Fix all and release only the lipid tails
model_sel.atomsel.set('occupancy', 1)
Selection('lipid and not name N "C1[1-5]" "H[1-5]\d" P1 "O."', model).atomsel.set('occupancy', 0)
model_sel.atomsel.write(FORMAT_PDB, LIPID_TAILS_PDB)

# Protein constraints
model_sel.atomsel.set('occupancy', 0)
Selection('protein and noh', model).atomsel.set('occupancy', 1)
model_sel.atomsel.write(FORMAT_PDB, PROTEIN_PDB)

# Protein CA constraints
model_sel.atomsel.set('occupancy', 0)
Selection('protein and name CA', model).atomsel.set('occupancy', 1)
model_sel.atomsel.write(FORMAT_PDB, BACKBONE_CA_PDB)

quit()
