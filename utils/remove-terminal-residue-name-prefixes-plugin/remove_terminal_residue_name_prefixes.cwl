#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Remove terminal residue name prefixes

doc: |-
  Remove terminal residue name prefixes

baseCommand: ['python', '-m', 'polus.mm.utils.remove_terminal_residue_name_prefixes']

hints:
  DockerRequirement:
    dockerPull: ndonyapour/remove_terminal_residue_name_prefixes

inputs:
  input_pdb_path:
    label: Path to the input file
    doc: |-
      Path to the input file
    type: File
    format: edam:format_1476 # pdb
    inputBinding:
      prefix: --input_pdb_path

  output_pdb_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: string
    format: edam:format_1476 # pdb
    inputBinding:
      prefix: --output_pdb_path
    default: system.pdb

outputs:
  output_pdb_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    format: edam:format_1476 # pdb
    outputBinding:
      glob: $(inputs.output_pdb_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
