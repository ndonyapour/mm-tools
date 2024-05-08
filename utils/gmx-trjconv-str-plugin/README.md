# gmx_trjconv_str (0.1.0)

Wrapper of the GROMACS trjconv module for converting between GROMACS compatible structure file formats and/or extracting a selection of atoms.

## Options

This plugin takes     5     input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_crd_path | Path to the input structure file, Type: string, File type: input, Accepted formats: xtc, trr, cpt, gro, g96, pdb, tng, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/trajectory.trr | Input | File | File |
| input_top_path | Path to the GROMACS input topology file, Type: string, File type: input, Accepted formats: tpr, gro, g96, pdb, brk, ent, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/topology.tpr | Input | File | File |
| output_str_path | Path to the output file, Type: string, File type: output, Accepted formats: pdb, xtc, trr, cpt, gro, g96, tng, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/reference/gromacs/ref_trjconv.str.pdb | Input | string | string |
| input_index_path | Path to the GROMACS index file, Type: string, File type: input, Accepted formats: ndx, Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/index.ndx | Input | File | File |
| config | Advanced configuration options for biobb_analysis GMXTrjConvStr. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_analysis GMXTrjConvStr documentation: https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_trjconv_str | Input | string | string |
| output_str_path | Path to the output file | Output | File | File |
