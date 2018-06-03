"""Build a model of a TRP channel."""
import tempfile

import VMD
from pyvmd import measure
from pyvmd.atoms import Selection
from pyvmd.molecules import Molecule, FORMAT_PDB, FORMAT_PSF

TEMPLATE_PDB = '../structures/3j5p.pdb'
TMP_DIR = tempfile.mkdtemp(prefix='vmd_')
TOPOLOGY = '/usr/lib/vmd/plugins/noarch/tcl/readcharmmtop1.1/top_all27_prot_lipid_na.inp'

# Load template
template = Molecule.create()
template.load(TEMPLATE_PDB)

# Center the template
center = measure.center(Selection('all', template))
Selection('all', template).atomsel.moveby((-center[0], -center[1], -center[2]))
# Rotate the template, so model takes less space
VMD.evaltcl('[atomselect %d all] move [transaxis z -35]' % template.molid)

# Save individual chains
for chain in 'ABCD':
    print "Saving monomer %s" % chain
    Selection('chain {} and resid 393 to 502'.format(chain), template).atomsel.write(FORMAT_PDB, '{}/{}1.pdb'.format(TMP_DIR, chain))
    Selection('chain {} and resid 508 to 719'.format(chain), template).atomsel.write(FORMAT_PDB, '{}/{}2.pdb'.format(TMP_DIR, chain))

# Load psfgen and topology
VMD.evaltcl("package require psfgen")
VMD.evaltcl("topology %s" % TOPOLOGY)
# Aliasing
VMD.evaltcl("pdbalias residue HIS HSE")
VMD.evaltcl("pdbalias atom ILE CD1 CD")
# Alias terminal oxygen for nitrogen from CT3 patch
VMD.evaltcl("pdbalias atom ALA OXT NT")

# Create segments
for chain in 'ABCD':
    VMD.evaltcl("segment %s1 {pdb %s/%s1.pdb; first ACE; last CT3}" % (chain, TMP_DIR, chain))
    VMD.evaltcl("coordpdb %s/%s1.pdb %s1" % (TMP_DIR, chain, chain))
    VMD.evaltcl("segment %s2 {pdb %s/%s2.pdb; first ACE; last CT3}" % (chain, TMP_DIR, chain))
    VMD.evaltcl("coordpdb %s/%s2.pdb %s2" % (TMP_DIR, chain, chain))

# Add 2 ions into selection filter
VMD.evaltcl('segment NA {residue 1 SOD; residue 2 SOD}')
VMD.evaltcl('coord NA 1 SOD {0 0 -33.5}')
VMD.evaltcl('coord NA 2 SOD {0 0 -37}')

# Guess coordinates of missing atoms
VMD.evaltcl("guesscoord")

# Save model with PSF
VMD.evaltcl("writepdb {}/model.pdb".format(TMP_DIR))
VMD.evaltcl("writepsf {}/model.psf".format(TMP_DIR))

model = Molecule.create()
model.load("{}/model.psf".format(TMP_DIR))
model.load("{}/model.pdb".format(TMP_DIR))

# Create bilayer
VMD.evaltcl('package require membrane')
VMD.evaltcl('membrane -l POPC -x 130 -y 130 -o %s/membrane -t c27' % TMP_DIR)

membrane = Molecule.create()
membrane.load('%s/membrane.psf' % TMP_DIR)
membrane.load('%s/membrane.pdb' % TMP_DIR)

# Move membrane to the center of the transmembrane part of the protein
trans_model_sel = Selection('resid 393 to 719', model).atomsel  # Transmembrane part
membrane_sel = Selection('all', membrane).atomsel
membrane_sel.moveby([c - m for m, c in zip(membrane_sel.center(), trans_model_sel.center())])

# Find conflicts between water
water_sel = Selection('water', membrane).atomsel
bad_atom_ids = set(water_sel.contacts(water_sel, 0.5)[0])
bad_str = ' '.join(str(i) for i in bad_atom_ids)
print "Found %d bad atoms (water)" % len(bad_atom_ids)
# Find water conflicting with lipids in membrane
bad_water = Selection('water and within 1.0 of not water', membrane)
bad_atom_ids.update({a.index for a in bad_water})
bad_str = ' '.join(str(i) for i in bad_atom_ids)
print "Found %d bad atoms (water-lipids)" % len(bad_atom_ids)
# Find confilicting lipids
model_sel = Selection('not hydrogen', model).atomsel
membrane_sel = Selection('not hydrogen and not same residue as (index %s)' % bad_str, membrane).atomsel
bad_atom_ids.update(membrane_sel.contacts(model_sel, 1.0)[0])
bad_str = ' '.join(str(i) for i in bad_atom_ids)
print "Found %d bad atoms (protein-lipids)" % len(bad_atom_ids)
# Find confilicting lipids - more strict, checks hydrogens, but in close range
model_sel = Selection('all', model).atomsel
membrane_sel = Selection('not same residue as (index %s)' % bad_str, membrane).atomsel
bad_atom_ids.update(membrane_sel.contacts(model_sel, 0.5)[0])
bad_str = ' '.join(str(i) for i in bad_atom_ids)
print "Found %d bad atoms (protein-lipids-H)" % len(bad_atom_ids)
# Find conflicting water (with protein) - more strict, checks hydrogens
model_sel = Selection('all', model).atomsel
membrane_sel = Selection('water and not same residue as (index %s)' % bad_str, membrane).atomsel
bad_atom_ids.update(membrane_sel.contacts(model_sel, 1.0)[0])
bad_str = ' '.join(str(i) for i in bad_atom_ids)
print "Found %d bad atoms (protein-water)" % len(bad_atom_ids)
# Write the membrane
good_membrane = Selection('not same residue as (index %s or sqr(x) + sqr(y) < sqr(20))' % bad_str, membrane)
good_membrane.atomsel.write(FORMAT_PSF, '%s/good_membrane.psf' % TMP_DIR)
good_membrane.atomsel.write(FORMAT_PDB, '%s/good_membrane.pdb' % TMP_DIR)

VMD.evaltcl('resetpsf')

# Merge protein and membrane
VMD.evaltcl('package require psfgen')
VMD.evaltcl("topology %s" % TOPOLOGY)

VMD.evaltcl('readpsf %s/model.psf' % TMP_DIR)
VMD.evaltcl('coordpdb %s/model.pdb' % TMP_DIR)
VMD.evaltcl('readpsf %s/good_membrane.psf' % TMP_DIR)
VMD.evaltcl('coordpdb %s/good_membrane.pdb' % TMP_DIR)

# Protein and membrane
VMD.evaltcl('writepsf %s/prot_memb.psf' % TMP_DIR)
VMD.evaltcl('writepdb %s/prot_memb.pdb' % TMP_DIR)

# Solvate system
# Use the X and Y dimension from membrane
mem_min, mem_max = Selection('water', membrane).atomsel.minmax()
# Compute Z dimension from protein and add 10 A water layer to top and bottom
prot_min, prot_max = Selection('protein', model).atomsel.minmax()
VMD.evaltcl('package require solvate')
# -minmax {{xmin ymin zmin} {xmax ymax zmax}}
data = {
    'tmp': TMP_DIR,
    'xmin': mem_min[0],
    'ymin': mem_min[1],
    'zmin': prot_min[2] - 15,
    'xmax': mem_max[0],
    'ymax': mem_max[1],
    'zmax': prot_max[2] + 15,
}
VMD.evaltcl('solvate %(tmp)s/prot_memb.psf %(tmp)s/prot_memb.pdb -o %(tmp)s/solvate -b 1.5 '
            '-minmax {{%(xmin)f %(ymin)f %(zmin)f} {%(xmax)f %(ymax)f %(zmax)f}}' % data)

# Remove extra water added by solvate, keep water from membrane
model = Molecule.create()
model.load('%s/solvate.psf' % TMP_DIR)
model.load('%s/solvate.pdb' % TMP_DIR)
low, high = Selection('lipid and element P', model).atomsel.minmax()
model_sel = Selection(
    'not same residue as (segid "WT.*" and z > %f and z < %f and not sqr(x)+sqr(y) <= sqr(10))' % (low[2], high[2]),
    model,
)
model_sel.atomsel.write(FORMAT_PSF, '%s/initial.psf' % TMP_DIR)
model_sel.atomsel.write(FORMAT_PDB, '%s/initial.pdb' % TMP_DIR)


# Neutralize system and add ions to make a 0.2 M NaCl
VMD.evaltcl('package require autoionize')
VMD.evaltcl('autoionize -psf %(tmp)s/initial.psf -pdb %(tmp)s/initial.pdb -o initial -sc 0.2' % {'tmp': TMP_DIR})

# Measurements
model = Molecule.create()
model.load('%s/initial.psf' % TMP_DIR)
model.load('%s/initial.pdb' % TMP_DIR)
# Measure only water as it is box-ish, lipids stick outside this box, minimization will take care of those.
box_sel = Selection('water', model)
print 'Center:', measure.center(box_sel)
# It is recommended to add 1 A in each dimension to get the box size
print 'Box size:', [h - l + 1 for l, h in zip(*box_sel.atomsel.minmax())]

quit()
