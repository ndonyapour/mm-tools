"""Package entrypoint for the {{cookiecutter.package_name}} package."""

# Base packages
import logging
from os import environ
import json

import typer
from {{cookiecutter.plugin_package}}.{{cookiecutter.package_name}} import (
    {{cookiecutter.package_name}},
)

logging.basicConfig(
    format="%(asctime)s - %(name)-8s - %(levelname)-8s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
POLUS_LOG = getattr(logging, environ.get("POLUS_LOG", "INFO"))
logger = logging.getLogger("{{cookiecutter.plugin_package}}.")
logger.setLevel(POLUS_LOG)

app = typer.Typer(help="{{cookiecutter.plugin_name}}.")

inputs = [
   "{{cookiecutter.inputs}}",
    "{{cookiecutter.outputs}}"
]

# convert each json string in the list to a dictionary via json.loads
inputs = [json.loads(i.replace("'", "\"")) for i in inputs]
# combine the dictionaries in the list to a single dictionary
inputs = {k: v for i in inputs for k, v in i.items()}

def get_options():
    """Returns a dictionary of named options."""
    options = {}
    for key, value in inputs.items():
        options[key] = typer.Option(..., help=value["description"])
    return options


@app.command()
def main(**kwargs):
    {{cookiecutter.plugin_description}}
    logger.info("Inputs:")
    for key, value in kwargs.items():
        value = inputs[key]
        logger.info(f"{key}: {value}")

    # Perform your logic with the dynamic inputs
    {{cookiecutter.package_name}}(**kwargs)


if __name__ == "__main__":
    options = get_options()

    for option_name, option in options.items():
        app.command(name=option_name)(option)

    main(**options)
