"""Read the inputs and outputs from a CWL file and insert them into the cookiecutter.json file."""
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
        cwl_content = yaml.load(file, Loader=yaml.FullLoader)

    parsed_data = {}
    for key in keys:
        parsed_data[key] = {}
        inputs = cwl_content.get(key, {})
        for input_name, input_data in inputs.items():
            input_type = input_data.get("type", "")
            if input_type == "string?":  # then dont include since its for edge cases
                continue
            parsed_data[key][input_name] = {"type": input_type}
            description = input_data.get("doc", "")
            if description:
                parsed_data[key][input_name]["description"] = description
            title = input_data.get("label", "")
            if title:
                parsed_data[key][input_name]["title"] = title

    return parsed_data


def transform_cwl(cwl: dict) -> dict:
    """Transform the cwl dictionary to include the python type and required fields.

    Args:
        cwl (Dict): The cwl dictionary.

    Returns:
        Dict: The transformed cwl dictionary.
    """
    for key in cwl:
        for input_data in cwl[key].values():
            if "title" not in input_data:
                input_data["title"] = ""
            if "description" not in input_data:
                input_data["description"] = ""
            if input_data["type"] == "File":
                input_data["python_type"] = "Path"
                input_data["type"] = "collection"
            else:
                input_data["python_type"] = ""
                input_data["type"] = ""
            input_data["required"] = "False"

    return cwl


def insert_inputs_outputs_cookiecutter(
    cookiecutter: dict,
    transformed: dict,
    cookiecutter_path: Path,
    keys: list[str],
) -> None:
    """Insert the inputs and outputs from the CWL file into the cookiecutter.json file.

    Args:
        cookiecutter (Dict): The cookiecutter dictionary.
        transformed (Dict): The transformed CWL dictionary.
        cookiecutter_path (Path): The path to the cookiecutter.json file.
        keys (List[str]): The keys to insert into the cookiecutter dictionary.
    """
    for key in keys:
        path_exists = any(
            input_data.get("python_type") == "Path"
            for input_data in transformed[key].values()
        )
        if not path_exists:
            # keep base template inpdir/outdir
            # only if there isnt any Path type
            for input_name, input_data in transformed[key].items():
                if input_name not in cookiecutter[key]:
                    cookiecutter[key][input_name] = input_data
        else:
            cookiecutter[key] = transformed[key]

    with cookiecutter_path.open("w", encoding="utf-8") as file:
        json.dump(cookiecutter, file, indent=4)


keys = ["inputs", "outputs"]
CWL_FILE = sys.argv[1]
result = parse_cwl(CWL_FILE, keys)
transformed_result = transform_cwl(result)
json_path = Path("cookiecutter.json")
with json_path.open(encoding="utf-8") as json_file:
    cookiecutter_data = json.load(json_file)
insert_inputs_outputs_cookiecutter(
    cookiecutter_data,
    transformed_result,
    json_path,
    keys,
)
