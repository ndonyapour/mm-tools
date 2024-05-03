#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

label: Returns a dictionary of the given arguments as a JSON-encoded string.
doc: |-
  Returns a dictionary of the given arguments as a JSON-encoded string.

baseCommand: echo # Anything, unused

requirements:
  InlineJavascriptRequirement: {}

inputs:
  water_type:
    label: Water molecule type. Values spc, spce, tip3p, tip4p, tip5p, tips3p.
    doc: |-
      Water molecule type. Values spc, spce, tip3p, tip4p, tip5p, tips3p.
    type: string
    format:
    - edam:format_2330

  forcefield:
    label: Force field to be used during the conversion. Values gromos45a3, charmm27, gromos53a6, amber96, amber99, gromos43a2, gromos54a7, gromos43a1, amberGS, gromos53a5, amber99sb, amber03, amber99sb-ildn.

    doc: |-
       Force field to be used during the conversion. Values gromos45a3, charmm27, gromos53a6, amber96, amber99, gromos43a2, gromos54a7, gromos43a1, amberGS, gromos53a5, amber99sb, amber03, amber99sb-ildn.
    type: string
    format: edam:format_2330

  ignh:
    label: Should pdb2gmx ignore the hidrogens in the original structure.
    doc: |-
       Should pdb2gmx ignore the hidrogens in the original structure.
    type: boolean
    format: edam:format_2330

  merge:
    label: Merge all chains into a single molecule.
    doc: |-
      Merge all chains into a single molecule.
    type: boolean
    format: edam:format_2330

# TODO: his

outputs:
  output_config_string:
    label: A dictionary of the given arguments as a JSON-encoded string.
    doc: |-
      A dictionary of the given arguments as a JSON-encoded string.
    type: string
    #format: edam:format_2330 # "'str' object does not support item assignment""
    outputBinding:
      outputEval: |
        ${
          var config = {};
          config["water_type"] = inputs.water_type;
          config["force_field"] = inputs.forcefield; // Note underscore
          config["ignh"] = inputs.ignh;
          config["merge"] = inputs.merge;
          return JSON.stringify(config);
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
