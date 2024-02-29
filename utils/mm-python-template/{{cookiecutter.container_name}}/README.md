# {{cookiecutter.plugin_name}} ({{cookiecutter.plugin_version}})

{{cookiecutter.plugin_description}}

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

