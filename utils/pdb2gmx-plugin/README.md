# pdb2gmx (0.1.0)

Wrapper class for the GROMACS pdb2gmx module.

## Options

This plugin takes 4 input arguments and 2 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_pdb_path | Path to the input PDB file, Type: string, File type: input, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs/egfr.pdb | Input | File | File |
| output_crd_path | Path to the output GRO file, Type: string, File type: output, Accepted formats: gro, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs/ref_pdb2gmx.gro | Input | string | string |
| output_top_zip_path | Path the output TOP topology in zip format, Type: string, File type: output, Accepted formats: zip, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs/ref_pdb2gmx.zip | Input | string | string |
| config | Advanced configuration options for biobb_gromacs Pdb2gmx. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_gromacs Pdb2gmx documentation: https://biobb-gromacs.readthedocs.io/en/latest/gromacs.html#module-gromacs.pdb2gmx | Input | string | string |
| output_crd_path | Path to the output GRO file | Output | File | File |
| output_top_zip_path | Path the output TOP topology in zip format | Output | File | File |
