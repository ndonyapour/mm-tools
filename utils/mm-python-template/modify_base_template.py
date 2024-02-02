"""Modify the base template to add input and output options dynamically."""
import json
from pathlib import Path

# Load cookiecutter.json
json_path = Path("cookiecutter.json")

with json_path.open() as json_file:
    cookiecutter_data = json.load(json_file)


def add_readme_options_dynamically(cookiecutter_json: dict) -> None:
    """Add the Options section to the README.md file dynamically.

    Args:
        cookiecutter_json (Dict): The cookiecutter.json file as a dictionary.
    """
    # Create the Options section dynamically
    options_section = "## Options\n\nThis plugin takes {} input arguments\
    and {} output argument:\n\n".format(
        len(cookiecutter_json["inputs"]),
        len(cookiecutter_json["outputs"]),
    )

    # Add the table headers
    options_section += (
        "| Name          | Description             | I/O    | Type   | Default |\n"
    )
    options_section += (
        "|---------------|-------------------------|--------|--------|---------|\n"
    )

    # Add rows for input arguments
    for input_name, input_data in cookiecutter_json["inputs"].items():
        options_section += "| {} | {} | {} | {} | {} |\n".format(
            input_name,
            input_data["description"],
            "Input",
            input_data["type"],
            input_data["type"],
        )

    # Add rows for output arguments
    for output_name, output_data in cookiecutter_json["outputs"].items():
        options_section += "| {} | {} | {} | {} | {} |\n".format(
            output_name,
            output_data["description"],
            "Output",
            output_data["type"],
            output_data["type"],
        )
    path = Path("{{cookiecutter.container_name}}") / "README.md"

    with path.open("a") as updated_readme_file:
        updated_readme_file.write(options_section)


def add_run_plugin_options_dynamically(cookiecutter_json: dict) -> None:
    """Add the input and output options to the run-plugin.sh file.

    Args:
        cookiecutter_json (Dict): The cookiecutter.json file as a dictionary.
    """
    docker_string = ""
    dict_items = ["inputs", "outputs"]
    for item in dict_items:
        for key in cookiecutter_json[item]:
            if key == "inpdir":
                value_string = "${container_input_dir}"
            elif key == "outdir":
                value_string = "${container_output_dir}"
            else:
                value_string = f"${{env_{key}}}"
            docker_string += f" --{key} {value_string}"

    path = Path("{{cookiecutter.container_name}}") / "run-plugin.sh"

    with path.open("a") as updated_run_plugin_file:
        updated_run_plugin_file.write(docker_string)


def add_function_arguments_dynamically(cookiecutter_json: dict) -> None:
    """Add the function arguments to the plugin file.

    Args:
        cookiecutter_json (Dict): The cookiecutter.json file as a dictionary.
    """
    package_name = cookiecutter_json["plugin_package"].split(".")[-1]
    def_string = f"def {package_name}("
    dict_items = ["inputs", "outputs"]
    for item in dict_items:
        for key, value in cookiecutter_json[item].items():
            def_string += f"{key} : {value['python_type']}, "
    def_string = def_string[:-2] + "):"
    # now add the docstring
    docstring = f"    '''{cookiecutter_json['plugin_name']}.\n\n    Args:\n"
    for item in dict_items:
        for key, value in cookiecutter_json[item].items():
            docstring += f"        {key}: {value['description']}\n"
    path = (
        Path("{{cookiecutter.container_name}}")
        / "src"
        / "{{cookiecutter.package_folders}}"
        / "{{cookiecutter.package_name}}.py"
    )

    with path.open("a") as file:
        file.write(
            def_string
            + "\n"
            + docstring
            + "    Returns:\n        None\n    '''\n    pass\n",
        )


add_readme_options_dynamically(cookiecutter_data)
add_run_plugin_options_dynamically(cookiecutter_data)
add_function_arguments_dynamically(cookiecutter_data)
