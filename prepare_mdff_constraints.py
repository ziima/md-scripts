"""Generate a constraints for MDFF simulation."""
import VMD
from pyvmd.molecules import Molecule

PSF = 'initial.psf'
PDB = 'initial.pdb'
MDFF_TEMPLATE_PSF = '../template.psf'
MDFF_TEMPLATE_PDB = '../template.pdb'
MDFF_MAP = 'mdff-map.situs'
MDFF_POTENTIAL = 'mdff-potential.dx'
MDFF_WEIGHTS = 'mdff-weights.pdb'
MDFF_CONSTRAINT_HBONDS = 'mdff-hbonds.dat'
MDFF_CONSTRAINT_CISPEPTIDE = 'mdff-cispeptide.dat'
MDFF_CONSTRAINT_CHIRALITY = 'mdff-chirality.dat'

model = Molecule.create()
model.load(PSF)
model.load(PDB)

template = Molecule.create()
template.load(MDFF_TEMPLATE_PSF)
template.load(MDFF_TEMPLATE_PDB)

# Load mdff
VMD.evaltcl('package require mdff')
# Generate simulated MDFF map
VMD.evaltcl('mdff sim [atomselect {molid} protein] -res 5 -o {outfile}'.format(molid=template.molid, outfile=MDFF_MAP))
# Convert map to potential
VMD.evaltcl('mdff griddx -i {infile} -o {outfile}'.format(infile=MDFF_MAP, outfile=MDFF_POTENTIAL))
# Generate scaling factors
VMD.evaltcl('mdff gridpdb -psf {psf} -pdb {pdb} -o {outfile}'.format(psf=PSF, pdb=PDB, outfile=MDFF_WEIGHTS))

# Generate additional constraints to prevent overfitting

# Prepare secondary structure restraints
VMD.evaltcl('package require ssrestraints')
VMD.evaltcl('ssrestraints -psf {psf} -pdb {pdb} -o {outfile} -hbonds'.format(psf=PSF, pdb=PDB,
                                                                             outfile=MDFF_CONSTRAINT_HBONDS))

# Prepare cis peptide restraints
VMD.evaltcl('package require cispeptide')
VMD.evaltcl('cispeptide restrain -mol {molid} -o {outfile}'.format(molid=model.molid,
                                                                   outfile=MDFF_CONSTRAINT_CISPEPTIDE))

# Prepare chirality restraints
VMD.evaltcl('package require chirality')
VMD.evaltcl('chirality restrain -mol {molid} -o {outfile}'.format(molid=model.molid, outfile=MDFF_CONSTRAINT_CHIRALITY))

quit()
