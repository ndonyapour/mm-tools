# editconf (0.1.0)

Wrapper class for the GROMACS editconf module.

## Options

This plugin takes 3 input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_crd_path | Path to the input GRO file, Type: string, File type: input, Accepted formats: gro, pdb, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs/editconf.gro | Input | File | File |
| output_crd_path | Path to the output GRO file, Type: string, File type: output, Accepted formats: pdb, gro, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs/ref_editconf.gro | Input | string | string |
| config | Advanced configuration options for biobb_gromacs Editconf. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_gromacs Editconf documentation: https://biobb-gromacs.readthedocs.io/en/latest/gromacs.html#module-gromacs.editconf | Input | string | string |
| output_crd_path | Path to the output GRO file | Output | File | File |
