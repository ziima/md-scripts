"""Compute distances between defined atoms pairs."""
import logging

from pyvmd.analyzer import Analyzer
from pyvmd.collectors import DistanceCollector
from pyvmd.datasets import DataSet
from pyvmd.molecules import Molecule

logging.basicConfig(level=logging.INFO)


PARAM_FILE = '../initial.psf'
PREFIXES = ['simulation_A', 'simulation_B']

COLLECTOR_1 = DistanceCollector('resname A2F and name N52B', 'protein and resid 1589 and name CD',
                                name="N52-GLU1589")
COLLECTOR_2 = DistanceCollector('resname A2F and name N49B', 'protein and resid 1534 and name CD',
                                name="N49-GLU1534")

MOLECULE = Molecule.create()
MOLECULE.load(PARAM_FILE)

for prefix in PREFIXES:
    trajectory = '../%s/%s.dcd' % (prefix, prefix)

    dset = DataSet()
    dset.add_collector(COLLECTOR_1)
    dset.add_collector(COLLECTOR_2)

    analyzer = Analyzer(MOLECULE, [trajectory])
    analyzer.add_dataset(dset)
    analyzer.analyze()

    dset.write('distances/%s.dat' % prefix)

quit()
