# config_tag_pdb2gmx (0.1.0)

Returns a dictionary of the given arguments as a JSON-encoded string.

## Options

This plugin takes     4     input arguments and 1 output argument:

| Name          | Description             | I/O    | Type   | Default |
|---------------|-------------------------|--------|--------|---------|
| water_type | Water molecule type. Values spc, spce, tip3p, tip4p, tip5p, tips3p. | Input | string | string |
| forcefield | Force field to be used during the conversion. Values gromos45a3, charmm27, gromos53a6, amber96, amber99, gromos43a2, gromos54a7, gromos43a1, amberGS, gromos53a5, amber99sb, amber03, amber99sb-ildn.  | Input | string | string |
| ignh | Should pdb2gmx ignore the hidrogens in the original structure. | Input | boolean | boolean |
| merge | Merge all chains into a single molecule. | Input | boolean | boolean |
| output_config_string | A dictionary of the given arguments as a JSON-encoded string. | Output | string | string |
