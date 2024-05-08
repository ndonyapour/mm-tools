# solvate (0.1.0)

Wrapper of the GROMACS solvate module.

## Options

This plugin takes     6     input arguments and 2 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_solute_crd_path | Path to the input GRO file, Type: string, File type: input, Accepted formats: gro, pdb, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs/solvate.gro | Input | File | File |
| output_crd_path | Path to the output GRO file, Type: string, File type: output, Accepted formats: gro, pdb, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs/ref_solvate.gro | Input | string | string |
| input_top_zip_path | Path the input TOP topology in zip format, Type: string, File type: input, Accepted formats: zip, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs/solvate.zip | Input | File | File |
| output_top_zip_path | Path the output topology in zip format, Type: string, File type: output, Accepted formats: zip, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs/ref_solvate.zip | Input | string | string |
| input_solvent_crd_path | (spc216.gro) Path to the GRO file containing the structure of the solvent, Type: string, File type: input, Accepted formats: gro, Example file: null | Input | File | File |
| config | Advanced configuration options for biobb_gromacs Solvate. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_gromacs Solvate documentation: https://biobb-gromacs.readthedocs.io/en/latest/gromacs.html#module-gromacs.solvate | Input | string | string |
| output_crd_path | Path to the output GRO file | Output | File | File |
| output_top_zip_path | Path the output topology in zip format | Output | File | File |
