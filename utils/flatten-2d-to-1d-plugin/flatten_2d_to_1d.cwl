#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: flatten a 2D array to 1D array
doc: |-
  flatten a 2D array to 1D array

baseCommand: python3

requirements:
  InlineJavascriptRequirement: {}


inputs:
  input_2d_array:
    label: Input 2D array # type:
    doc: |-
      Input 2D array
      Type: string[][]
      File type: input
      Accepted formats: list[list[string]]
    type: ["null", {"type": "array", "items": {"type": "array", "items": "string"}}]
    format: edam:format_2330

outputs:
  output_1d_array:
    label: Output 1D array
    doc: |-
      Output 1D array
    type: ["null", {"type": "array", "items": "string"}]
    outputBinding:
      outputEval: |
        ${
          // if (inputs.input_2D_array == null) {
          //   return null;
          // }
          var pdbids_1d = [];
          for (var i = 0; i < inputs.input_2d_array.length; i++) {
            var pdbids = inputs.input_2d_array[i];
             if (pdbids != null) { // Check if the 1D array is not null
              for (var j = 0; j < pdbids.length; j++){
                if (pdbids[j] != null){
                  pdbids_1d.push(pdbids[j]);
                }
              }
             }
          }
        return pdbids_1d;
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
