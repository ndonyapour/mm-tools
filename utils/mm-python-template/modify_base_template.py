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
    # add top level imports
    import_string = "from typing import List\nfrom pathlib import Path\n"
    # create definition string
    for item in dict_items:
        for key, value in cookiecutter_json[item].items():
            # use filepattern in __main__ to pass multiple files form inpdir
            if key == "inpdir":
                continue
            if "python_type" in value:
                def_string += f"{key} : {value['python_type']}, "
            else:
                if "_pattern" in key:
                    original_key = key.replace("_pattern", "")
                    new_python_type = "List[Path]"
                    def_string += f"{original_key}: {new_python_type}, "
    def_string = def_string[:-2] + "):"
    # now add the docstring
    docstring = f"    '''{cookiecutter_json['plugin_name']}.\n\n    Args:\n"
    for item in dict_items:
        for key, value in cookiecutter_json[item].items():
            if key == "inpdir":
                continue
            if "python_type" in value:
                docstring += f"        {key}: {value['description']}\n"
            else:
                if "_pattern" in key:
                    original_key = key.replace("_pattern", "")
                    docstring += f"        {original_key}: {value['description']}\n"
    path = (
        Path("{{cookiecutter.container_name}}")
        / "src"
        / "{{cookiecutter.package_folders}}"
        / "{{cookiecutter.package_name}}.py"
    )

    with path.open("a") as file:
        file.write(
            import_string
            + "\n"
            + def_string
            + "\n"
            + docstring
            + "    Returns:\n        None\n    '''\n    pass\n",
        )


def add_main_function_dynamically(cookiecutter_json: dict) -> None:
    """Add the main function to the __main__.py file.

    Args:
        cookiecutter_json (Dict): The cookiecutter.json file as a dictionary.
    """
    # Create the main function string
    main_function_string = "\n\n@app.command()\ndef main("

    dict_items = ["inputs", "outputs"]
    for item in dict_items:
        for key, value in cookiecutter_json[item].items():
            main_function_string += (
                f"\n    {key}: {value['python_type']} = typer.Option("
            )
            main_function_string += "\n        ...,"
            main_function_string += f"\n        '--{key}',"
            main_function_string += f"\n        help='{value['description']}',"
            main_function_string += "\n    ),"

    # Add the main function body
    main_function_string += "\n) -> None:"
    main_function_string += f'\n    """{cookiecutter_json["plugin_name"]}."""'
    main_function_string += '\n    logger.info(f"pattern: {input_pdb_path_pattern}")'
    for item in dict_items:
        for key, value in cookiecutter_json[item].items():
            main_function_string += f'\n    logger.info(f"{key}: {{{key}}}")'

    main_function_string += "\n    file_paths = []"
    main_function_string += (
        "\n    files = fp.FilePattern(inpdir, input_pdb_path_pattern)"
    )
    main_function_string += "\n    for file in files():"
    main_function_string += "\n        file_path = file[-1][0]"
    main_function_string += "\n        file_paths.append(file_path)"
    main_function_string += "\n    {{cookiecutter.package_name}}(file_paths, outdir)\n"

    # Check if the __main__.py file exists
    main_file_path = (
        Path("{{cookiecutter.container_name}}")
        / "src"
        / "{{cookiecutter.package_folders}}"
        / "__main__.py"
    )

    with main_file_path.open("a") as file:
        file.write(main_function_string)

def add_ict_inputs_outputs_dynamically(cookiecutter_json: dict) -> None:
    """Add the input and output options to the ict.yaml file.

    Args:
        cookiecutter_json (Dict): The cookiecutter.json file as a dictionary.
    """
    path = Path("{{cookiecutter.container_name}}") / "ict.yaml"

    array = ["inputs", "outputs"]
    with path.open("a") as updated_ict_file:
        for item in array:
            updated_ict_file.write(f"{item}:\n")
            for key, value in cookiecutter_json[item].items():
                updated_ict_file.write(f"  - name: {key}\n")
                updated_ict_file.write(f"    required: true\n")
                updated_ict_file.write(f"    description: {value['description']}\n")
                updated_ict_file.write(f"    type: {value['type']}\n")

# need to handle cases where base_command is just python3
# calling some script already in container for example
ADD_SOURCE = True
base_command = cookiecutter_data["base_command"]
if len(base_command.split()) == 1:
    ADD_SOURCE = False

add_readme_options_dynamically(cookiecutter_data)
add_run_plugin_options_dynamically(cookiecutter_data)
add_ict_inputs_outputs_dynamically(cookiecutter_data)
if ADD_SOURCE:
    add_function_arguments_dynamically(cookiecutter_data)
    add_main_function_dynamically(cookiecutter_data)
    cookiecutter_data[
        "base_command"
    ] = f"python3 -m {cookiecutter_data['plugin_package']}"
    with json_path.open("w") as json_file:
        json.dump(cookiecutter_data, json_file, indent=4)

else:
    # in this case we actually dont want a src folder
    # need a way to tell cookiecutter to exclude src folder

    # In this case since calling some internal tool dont add additional
    # function arguments such as inpdir and outdir
    del cookiecutter_data["inputs"]["inpdir"]
    del cookiecutter_data["outputs"]["outdir"]
    with json_path.open("w") as json_file:
        json.dump(cookiecutter_data, json_file, indent=4)
