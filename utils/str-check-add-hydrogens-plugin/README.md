# str_check_add_hydrogens (0.1.0)

This class is a wrapper of the Structure Checking tool to add hydrogens to a 3D structure.

## Options

This plugin takes 3 input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_structure_path | Input structure file path, Type: string, File type: input, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/data/utils/str_no_H.pdb | Input | File | File |
| output_structure_path | Output structure file path, Type: string, File type: output, Accepted formats: pdb, pdbqt, Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/reference/utils/ref_str_H.pdbqt | Input | string | string |
| config | Advanced configuration options for biobb_structure_utils StrCheckAddHydrogens. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_structure_utils StrCheckAddHydrogens documentation: https://biobb-structure-utils.readthedocs.io/en/latest/utils.html#utils-str-check-add-hydrogens-module | Input | string | string |
| output_structure_path | Output structure file path | Output | File | File |
