"""Tests for remove_terminal_residue_name_prefixes."""
import sys
from pathlib import Path

from polus.mm.utils.remove_terminal_residue_name_prefixes.remove_terminal_residue_name_prefixes import (  # noqa: E501
    remove_terminal_residue_name_prefixes,
)

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_remove_terminal_residue_name_prefixes() -> None:
    """Test remove_terminal_residue_name_prefixes."""
    input_pdb_path = "protein_1e3g.pdb"
    input_pdb_path = str(Path(__file__).resolve().parent / Path(input_pdb_path))

    output_pdb_path = "protein.pdb"
    output_pdb_path = str(Path(__file__).resolve().parent / Path(output_pdb_path))

    remove_terminal_residue_name_prefixes(input_pdb_path, output_pdb_path)
    assert Path(output_pdb_path).exists()


def test_cwl_remove_terminal_residue_name_prefixes() -> None:
    """Test the remove_terminal_residue_name_prefixes function in cwltool."""
    cwl_file = Path("remove_terminal_residue_name_prefixes.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)
    input_pdb_path = "protein_1e3g.pdb"
    input_pdb_path = str(Path(__file__).resolve().parent / Path(input_pdb_path))
    input_to_props["input_pdb_path"]["path"] = input_pdb_path

    input_to_props["output_pdb_path"] = "output.pdb"
    input_yaml_path = Path("remove_terminal_residue_name_prefixes.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)
    assert Path("output.pdb").exists()
