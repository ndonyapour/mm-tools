"""Tests for str_check_add_hydrogens."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_str_check_add_hydrogens() -> None:
    """Test str_check_add_hydrogens."""
    cwl_file = Path("str_check_add_hydrogens.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)
    input_pdb_path = "protein_1e3g.pdb"
    input_pdb_path = str(Path(__file__).resolve().parent / Path(input_pdb_path))
    input_to_props["input_structure_path"]["path"] = input_pdb_path
    input_to_props[
        "config"
    ] = '{"charges": false, "mode": "auto", "ph": 7.4, "list": "", \
    "keep_canonical_resnames": false}'

    input_to_props["output_structure_path"] = "output.pdb"
    input_yaml_path = Path("str_check_add_hydrogens.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)

    assert Path("output.pdb").exists()
