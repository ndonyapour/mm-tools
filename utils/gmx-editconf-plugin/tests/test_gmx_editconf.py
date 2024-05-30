"""Tests for gmx_editconf."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_gmx_editconf() -> None:
    """Test gmx_editconf."""
    cwl_file = Path("gmx_editconf.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)
    input_pdb_path = "1e3g.pdb"
    input_pdb_path = str(Path(__file__).resolve().parent / Path(input_pdb_path))
    input_to_props["input_crd_path"]["path"] = input_pdb_path
    input_to_props["align_principal_axes"] = 0
    input_to_props["box_type"] = "cubic"
    input_to_props["distance_to_molecule"] = 1.2
    input_to_props["output_crd_path"] = "system.g96"
    input_to_props["box_vector_lengths"] = [8.1, 8.1, 8.1]
    input_to_props["box_vector_angles"] = [90, 90, 90]

    input_yaml_path = Path("gmx_editconf.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)
    assert Path("system.g96").exists()
