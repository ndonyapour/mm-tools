"""Tests for fix_amides."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_fix_amides() -> None:
    """Test fix_amides."""
    cwl_file = Path("fix_amides.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)

    input_pdb_path = "4xk9.pdb"
    input_pdb_path = str(Path(__file__).resolve().parent / Path(input_pdb_path))

    input_to_props["input_pdb_path"]["path"] = input_pdb_path
    input_to_props["output_pdb_path"] = "protein_fix_amides.pdb"

    input_yaml_path = Path("fix_amides.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)

    assert Path("protein_fix_amides.pdb").exists()
