# convert_mol2 (0.1.0)

This class is a wrapper of the Open Babel tool.

## Reading inputs/outputs from .cwl files
This adds inputs/outputs from .cwl files into cookiecutter.json
`python read_cwl_inputs_outputs.py path_to_cwl_file.cwl`

## Modifying template files
To dynamically add inputs/outputs from cookiecutter.json to README.MD, __main__.py and plugin_package function
`python modify_base_template.py`

## Building

To build the Docker image for the conversion plugin, run `./build-docker.sh`.

## Install WIPP Plugin

If WIPP is running, navigate to the plugins page and add a new plugin. Paste the
contents of `plugin.json` into the pop-up window and submit.

## Options

This plugin takes 6 input arguments    and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| inpdir | Input file collection to be processed. | Input | path | path |
| first_molecule | Input Index of the first molecule (1-based), Type: string? | Input | string | string |
| last_molecule | Input Index of the last molecule (1-based), Type: string? | Input | string | string |
| output_mol2_path | Path to the output file, Type: string, File type: output, Accepted formats: mol2 | Input | string | string |
| arg1 | Additional arguments, Type: string? | Input | string | string |
| input_path_pattern | Path to the input file, Type: string, File type: input, Accepted formats: dat, ent, fa, fasta, gro, inp, log, mcif, mdl, mmcif, mol, mol2, pdb, pdbqt, png, sdf, smi, smiles, txt, xml, xtc, Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/data/babel/babel.smi | Input | string | string |
| outdir | Output collection. | Output | path | path |
