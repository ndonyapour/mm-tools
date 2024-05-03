"""Tests for config_tag_pdb2gmx."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_config_tag_pdb2gmx() -> None:
    """Test config_tag_pdb2gmx."""
    cwl_file = Path("config_tag_pdb2gmx.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)

    input_to_props["water_type"] = "spec"
    input_to_props["forcefield"] = "mber99sb-ildn"
    input_to_props["ignh"] = True
    input_to_props["merge"] = False
    input_yaml_path = Path("config_tag_pdb2gmx.yml")

    create_input_yaml(input_to_props, input_yaml_path)
    stdout, stderr = call_cwltool(cwl_file, input_yaml_path)
    assert "success" in stderr
