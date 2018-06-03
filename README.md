# md-scripts
Repository with various scripts used for molecular dynamics.
The scripts are usually designed with a single purpose and need to be modified for any other use.

Most of the scripts require [pyvmd](//github.com/ziima/pyvmd) library and they can be well combined with the [md-project-template](//github.com/ziima/md-project-template).

## List of scripts ##
### [download_trp_tails.py](download_trp_tails.py) ###
A Python script to download a sequences of human TRP channel tails from UniProt database.

### [build_trp_model.py](build_trp_model.py) ###
A Python VMD script to generate a model of a TRP channel from a PDB file.

### [prepare_trp_constraints.py](prepare_trp_constraints.py) ###
A Python VMD script to generate constraints for a simulation of a TRP channel.

### [prepare_mdff_constraints.py](prepare_mdff_constraints.py) ###
A Python VMD script to generate constraints for a MDFF simulation.

### [collect_distances.py](collect_distances.py) ###
A Python VMD script to collect distances of pairs of atoms over a trajectory.

### [parse_fep.py](parse_fep.py) ###
A Python VMD script to parse FEP simulation output.
