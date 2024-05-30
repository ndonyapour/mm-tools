"""Tests for gmx_trjconv_str."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_gmx_trjconv_str() -> None:
    """Test gmx_trjconv_str."""
    cwl_file = Path("gmx_trjconv_str.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)

    input_crd_path = "protein_1e3g.gro"
    input_crd_path = str(Path(__file__).resolve().parent / Path(input_crd_path))
    input_to_props["input_crd_path"]["path"] = input_crd_path

    input_to_props["input_top_path"]["path"] = input_crd_path

    input_to_props["output_str_path"] = "genion.pdb"
    del input_to_props["input_index_path"]

    input_yaml_path = Path("gmx_trjconv_str.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)
    assert Path("genion.pdb").exists()
