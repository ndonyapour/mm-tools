# fix_amides (0.1.0)

Fix amide groups from residues.

## Options

This plugin takes 3 input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_pdb_path | Input PDB file path, Type: string, File type: input, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_model/raw/master/biobb_model/test/data/model/5s2z.pdb | Input | File | File |
| output_pdb_path | Output PDB file path, Type: string, File type: output, Accepted formats: pdb, Example file: https://raw.githubusercontent.com/bioexcel/biobb_model/master/biobb_model/test/reference/model/output_amide_pdb_path.pdb | Input | string | string |
| config | Advanced configuration options for biobb_model FixAmides. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_model FixAmides documentation: https://biobb-model.readthedocs.io/en/latest/model.html#module-model.fix_amides | Input | string | string |
| output_pdb_path | Output PDB file path | Output | File | File |
