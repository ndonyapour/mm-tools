"""Package entrypoint for the remove_terminal_residue_name_prefixes package."""

# Base packages
import logging
from os import environ

import typer
from polus.mm.utils.remove_terminal_residue_name_prefixes.remove_terminal_residue_name_prefixes import (  # noqa: E501
    remove_terminal_residue_name_prefixes,
)

logging.basicConfig(
    format="%(asctime)s - %(name)-8s - %(levelname)-8s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
POLUS_LOG = getattr(logging, environ.get("POLUS_LOG", "INFO"))
logger = logging.getLogger("polus.mm.utils.remove_terminal_residue_name_prefixes.")
logger.setLevel(POLUS_LOG)

app = typer.Typer(help="remove_terminal_residue_name_prefixes.")


@app.command()
def main(
    input_pdb_path: str = typer.Option(
        ...,
        "--input_pdb_path",
        help="",
    ),
    output_pdb_path: str = typer.Option(
        ...,
        "--output_pdb_path",
        help="Path to the output file",
    ),
) -> None:
    """remove_terminal_residue_name_prefixes."""
    logger.info(f"input_pdb_path: {input_pdb_path}")
    logger.info(f"output_pdb_path: {output_pdb_path}")

    remove_terminal_residue_name_prefixes(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_pdb_path,
    )


if __name__ == "__main__":
    app()
