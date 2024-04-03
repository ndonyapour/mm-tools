"""Package entrypoint for the extract_protein_openmm package."""

# Base packages
import logging
from os import environ
from pathlib import Path

import typer
from polus.mm.utils.extract_protein_openmm.extract_protein_openmm import (
    extract_protein_openmm,
)

logging.basicConfig(
    format="%(asctime)s - %(name)-8s - %(levelname)-8s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
POLUS_LOG = getattr(logging, environ.get("POLUS_LOG", "INFO"))
logger = logging.getLogger("polus.mm.utils.extract_protein_openmm.")
logger.setLevel(POLUS_LOG)

app = typer.Typer(help="extract_protein_openmm.")



@app.command()
def main(
    input_pdb_path: str = typer.Option(
        ...,
        '--input_pdb_path',
        help='Input pdb file path, Type: string, File type: input, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/data/utils/cat_protein.pdb',
    ),
    output_pdb_path: str = typer.Option(
        ...,
        '--output_pdb_path',
        help='Output pdb file path, Type: string, File type: output, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/reference/utils/ref_cat_pdb.pdb',
    ),
) -> None:
    """extract_protein_openmm."""

    output_pdb_path = str(output_pdb_path)
    logger.info(f"input_pdb_path: {input_pdb_path}")
    logger.info(f"output_pdb_path: {output_pdb_path}")

    extract_protein_openmm(input_pdb_path, output_pdb_path)


if __name__ == "__main__":
    app()