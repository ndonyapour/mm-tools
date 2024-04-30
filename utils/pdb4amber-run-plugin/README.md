# pdb4amber_run (0.1.0)

Wrapper of the AmberTools (AMBER MD Package) pdb4amber tool module.

## Options

This plugin takes 3  input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_pdb_path | Input 3D structure PDB file, Type: string, File type: input, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/pdb4amber/1aki_fixed.pdb | Input | File | File |
| output_pdb_path | Output 3D structure PDB file, Type: string, File type: output, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/reference/pdb4amber/structure.pdb4amber.pdb | Input | string | string |
| config | Advanced configuration options for biobb_amber.pdb4amber.pdb4amber_run Pdb4amberRun. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_amber.pdb4amber.pdb4amber_run Pdb4amberRun documentation: https://biobb-amber.readthedocs.io/en/latest/pdb4amber.html#module-pdb4amber.pdb4amber_run | Input | string | string |
| output_pdb_path | Output 3D structure PDB file | Output | File | File |
