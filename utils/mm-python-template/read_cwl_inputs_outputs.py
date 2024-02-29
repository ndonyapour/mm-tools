"""Read the inputs and outputs from a CWL file and add to cookiecutter.json."""
import json
import sys
from pathlib import Path

from ruamel.yaml import YAML


def parse_cwl(file_path: str, keys: list[str]) -> dict:
    """Parse the inputs and outputs from a CWL file.

    Args:
        file_path (str): The path to the CWL file.
        keys (List[str]): The keys to parse from the CWL file.

    Returns:
        Dict: A dictionary with the inputs and outputs.
    """
    yaml = YAML(typ="safe", pure=True)
    with Path(file_path).open(encoding="utf-8") as file:
        cwl_content = yaml.load(file)

    # string? used for explicit edges
    black_list_types = ["string?"]

    parsed_data = {}
    for key in keys:
        parsed_data[key] = {}
        inputs = cwl_content.get(key, {})
        # just keep all outputs in outdir
        # dont specify each output
        if key == "outputs":
            continue
        for input_name, input_data in inputs.items():
            input_type = input_data.get("type", "")
            if input_type in black_list_types:
                continue
            # want to autogenerate output filenames in unique fashion
            # that depends on the input filename so when scattering (not with cwltool)
            # there is no name collision
            if "output" in input_name and input_type == "string":
                continue
            parsed_data[key][input_name] = {"type": input_type}
            description = input_data.get("doc", "")
            if description:
                description = description.replace("\n", ", ")
                parsed_data[key][input_name]["description"] = description
            title = input_data.get("label", "")
            if title:
                parsed_data[key][input_name]["title"] = title
            default = input_data.get("default", "")
            if default:
                parsed_data[key][input_name]["default"] = default
            input_binding = input_data.get("inputBinding", {})
            if input_binding:
                prefix = input_binding.get("prefix", "")
                if prefix:
                    parsed_data[key][input_name]["prefix"] = prefix
                else:
                    parsed_data[key][input_name]["prefix"] = ""

    # determine plugin_name
    hints = cwl_content.get("hints", {})
    docker_requirement = hints.get("DockerRequirement", {})
    docker_pull = docker_requirement.get("dockerPull", "")
    if docker_pull:
        plugin_name = docker_pull.split("/")[-1]
        parsed_data["plugin_name"] = plugin_name
    else:
        docker_image_id = docker_requirement.get("dockerImageId", "")
        parsed_data["plugin_name"] = docker_image_id

    # determine baseCommand
    base_command = cwl_content.get("baseCommand", "")
    if base_command:
        if isinstance(base_command, list):
            base_command = " ".join(base_command)
        parsed_data["base_command"] = base_command

    return parsed_data


def transform_cwl(cwl: dict, keys: list[str]) -> dict:
    """Transform the cwl dictionary to include the python type and required fields.

    Args:
        cwl (Dict): The cwl dictionary.
        keys (List[str]): The keys to transform in the cwl dictionary.

    Returns:
        Dict: The transformed cwl dictionary.
    """
    python_transforms = {"string": "str", "boolean": "bool"}
    transformed_input_keys = []
    for key in cwl:
        if key in keys:
            for input_key, input_data in cwl[key].items():
                if "prefix" not in input_data:
                    input_data["prefix"] = f'--{input_key}'
                if "title" not in input_data:
                    input_data["title"] = ""
                if "description" not in input_data:
                    input_data["description"] = ""
                # if default exists, then required is False else True
                if "default" not in input_data:
                    input_data["required"] = "True"
                else:
                    input_data["required"] = "False"
                if input_data["type"] == "File":
                    # want to use filepattern for these inputs
                    input_data["type"] = "string"
                    input_data["python_type"] = "str"
                    if key == "inputs":
                        transformed_input_keys.append(input_key)
                else:
                    # if ? in type, then remove the ?
                    if "?" in input_data["type"]:
                        input_data["type"] = input_data["type"].replace("?", "")
                    input_data["python_type"] = input_data["type"]
                    if input_data["python_type"] in python_transforms:
                        input_data["python_type"] = python_transforms[
                            input_data["python_type"]
                        ]


    # now for all keys in transformed_input_keys, append _pattern to the name
    # want to use filepattern for these inputs
    for key in transformed_input_keys:
        cwl["inputs"][f"{key}_pattern"] = cwl["inputs"].pop(key)

    return cwl


def insert_inputs_outputs_cookiecutter(
    cookiecutter: dict,
    transformed: dict,
    cookiecutter_path: Path,
    keys: list[str],
) -> None:
    """Insert the inputs and outputs from the CWL file in cookiecutter.json.

    Args:
        cookiecutter (Dict): The cookiecutter dictionary.
        transformed (Dict): The transformed CWL dictionary.
        cookiecutter_path (Path): The path to the cookiecutter.json file.
        keys (List[str]): The keys to insert into the cookiecutter dictionary.
    """
    # add base_command to cookiecutter
    cookiecutter["base_command"] = transformed["base_command"]

    for key in keys:
        for input_name, input_data in transformed[key].items():
            if input_name not in cookiecutter[key]:
                cookiecutter[key][input_name] = input_data

    cookiecutter["plugin_name"] = transformed["plugin_name"]
    plugin_package = cookiecutter["plugin_package"]
    plugin_package = plugin_package.split(".")
    plugin_package.pop()
    plugin_package = ".".join(plugin_package)
    plugin_package = f"{plugin_package}.{transformed['plugin_name']}"
    cookiecutter["plugin_package"] = plugin_package

    with cookiecutter_path.open("w", encoding="utf-8") as file:
        json.dump(cookiecutter, file, indent=4)


keys = ["inputs", "outputs"]
CWL_FILE = sys.argv[1]
result = parse_cwl(CWL_FILE, keys)
transformed_result = transform_cwl(result, keys)
json_path = Path("cookiecutter.json")
with json_path.open(encoding="utf-8") as json_file:
    cookiecutter_data = json.load(json_file)
insert_inputs_outputs_cookiecutter(
    cookiecutter_data,
    transformed_result,
    json_path,
    keys,
)
