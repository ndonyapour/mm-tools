"""Tests for convert_xyz."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_convert_xyz() -> None:
    """Test convert_xyz."""
    cwl_file = Path("convert_xyz.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)
    input_gro_path = "protein_1e3g.gro"
    input_gro_path = str(Path(__file__).resolve().parent / Path(input_gro_path))
    input_to_props["input_path"]["path"] = input_gro_path

    input_to_props["output_xyz_path"] = "protein.xyz"
    input_yaml_path = Path("convert_xyz.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)

    assert Path("protein.xyz").exists()
