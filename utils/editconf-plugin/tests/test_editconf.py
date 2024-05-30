"""Tests for editconf."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_editconf() -> None:
    """Test editconf."""
    cwl_file = Path("editconf.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)
    input_pdb_path = "1e3g.pdb"
    input_pdb_path = str(Path(__file__).resolve().parent / Path(input_pdb_path))
    input_to_props["input_crd_path"]["path"] = input_pdb_path
    input_to_props["config"] = '{"box_type": "cubic", "distance_to_molecule": 1.2}'
    input_to_props["output_crd_path"] = "system.g96"

    input_yaml_path = Path("editconf.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)

    assert Path("system.g96").exists()
