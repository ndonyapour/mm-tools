# gmx_editconf (0.1.0)

Wrapper class for the GROMACS editconf module.

## Options

This plugin takes     7     input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| input_crd_path | Path to the input GRO file, Type: string, File type: input, Accepted formats: gro, pdb, Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/editconf.gro | Input | File | File |
| output_crd_path | Path to the output GRO file, Type: string, File type: output, Accepted formats: pdb, gro, Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_editconf.gro | Input | string | string |
| distance_to_molecule | Distance between the solute and the box | Input | float | float |
| box_type | Box type such as triclinic, cubic, dodecahedron, octahedron | Input | string | string |
| box_vector_lengths | Box vector lengths (a,b,c) | Input | ['null', {'type': 'array', 'items': 'float'}] | ['null', {'type': 'array', 'items': 'float'}] |
| box_vector_angles | Angles between the box vectors (bc,ac,ab) | Input | ['null', {'type': 'array', 'items': 'float'}] | ['null', {'type': 'array', 'items': 'float'}] |
| align_principal_axes |  | Input | int | int |
| output_crd_path | Path to the output GRO file | Output | File | File |
