"""Tests for genion."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_genion() -> None:
    """Test genion."""
    cwl_file = Path("genion.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)

    input_tpr_path = "1e3g_system.tpr"
    input_tpr_path = str(Path(__file__).resolve().parent / Path(input_tpr_path))
    input_to_props["input_tpr_path"]["path"] = input_tpr_path

    input_top_zip_path = "1e3g_system.zip"
    input_top_zip_path = str(Path(__file__).resolve().parent / Path(input_top_zip_path))
    input_to_props["input_top_zip_path"]["path"] = input_top_zip_path
    input_to_props["output_crd_path"] = "1e3g_genion.gro"
    del input_to_props["input_ndx_path"]

    input_yaml_path = Path("genion.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)
    assert Path("1e3g_genion.gro").exists()
