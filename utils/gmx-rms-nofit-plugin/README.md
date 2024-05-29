# gmx_rms_nofit (0.1.0)

Wrapper of the GROMACS rms module for performing a Root Mean Square deviation (RMSd) analysis from a given GROMACS compatible trajectory.

## Options

This plugin takes     4     input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_structure_path | Path to the input structure file, Type: string, File type: input, Accepted formats: tpr, gro, g96, pdb, brk, ent, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/topology.tpr | Input | File | File |
| input_traj_path | Path to the GROMACS trajectory file, Type: string, File type: input, Accepted formats: xtc, trr, cpt, gro, g96, pdb, tng, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/trajectory.trr | Input | File | File |
| output_xvg_path | Path to the XVG output file, Type: string, File type: output, Accepted formats: xvg, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/reference/gromacs/ref_rms.xvg | Input | string | string |
| input_index_path | Path to the GROMACS index file, Type: string, File type: input, Accepted formats: ndx, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/index.ndx | Input | File | File |
| output_xvg_path | Path to the XVG output file | Output | File | File |
