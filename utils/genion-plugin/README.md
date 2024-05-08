# genion (0.1.0)

Wrapper class for the GROMACS genion module.

## Options

This plugin takes     6     input arguments and 2 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_tpr_path | Path to the input portable run input TPR file, Type: string, File type: input, Accepted formats: tpr, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs/genion.tpr | Input | File | File |
| output_crd_path | Path to the input structure GRO file, Type: string, File type: output, Accepted formats: gro, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs/ref_genion.gro | Input | string | string |
| input_top_zip_path | Path the input TOP topology in zip format, Type: string, File type: input, Accepted formats: zip, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs/genion.zip | Input | File | File |
| output_top_zip_path | Path the output topology TOP and ITP files zipball, Type: string, File type: output, Accepted formats: zip, Example file: https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs/ref_genion.zip | Input | string | string |
| input_ndx_path | Path to the input index NDX file, Type: string, File type: input, Accepted formats: ndx, Example file: null | Input | File | File |
| config | Advanced configuration options for biobb_gromacs Genion. This should be passed as a string containing a dict. The possible options to include here are listed under 'properties' in the biobb_gromacs Genion documentation: https://biobb-gromacs.readthedocs.io/en/latest/gromacs.html#module-gromacs.genion | Input | string | string |
| output_crd_path | Path to the input structure GRO file | Output | File | File |
| output_top_zip_path | Path the output topology TOP and ITP files zipball | Output | File | File |
