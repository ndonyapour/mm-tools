# box (0.1.0)

This class sets the center and the size of a rectangular parallelepiped box around a set of residues or a pocket.

## Options

This plugin takes 3 input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_pdb_path | PDB file containing a selection of residue numbers or PQR file containing the pocket, Type: string, File type: input, Accepted formats: pdb, pqr, Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/utils/input_box.pqr | Input | File | File |
| output_pdb_path | PDB including the annotation of the box center and size as REMARKs, Type: string, File type: output, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/utils/ref_output_box.pdb | Input | string | string |
| config | Advanced configuration options for biobb_vs Box. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_vs Box documentation: https://biobb-vs.readthedocs.io/en/latest/utils.html#module-utils.box | Input | string | string |
| output_pdb_path | PDB including the annotation of the box center and size as REMARKs | Output | File | File |
