"""Tests for pymol_align_protein_ca."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_pymol_align_protein_ca() -> None:
    """Test pymol_align_protein_ca."""
    cwl_file = Path("pymol_align_protein_ca.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)

    input_1_path = "receptor.xyz"
    input_1_path = str(Path(__file__).resolve().parent / Path(input_1_path))
    input_to_props["input_1_path"]["path"] = input_1_path

    input_2_path = "pose_ligand.xyz"
    input_2_path = str(Path(__file__).resolve().parent / Path(input_2_path))
    input_to_props["input_2_path"]["path"] = input_2_path

    input_3_path = "npt.gro"
    input_3_path = str(Path(__file__).resolve().parent / Path(input_3_path))
    input_to_props["input_3_path"]["path"] = input_3_path

    input_4_path = "prod.trr"
    input_4_path = str(Path(__file__).resolve().parent / Path(input_4_path))
    input_to_props["input_4_path"]["path"] = input_4_path

    input_to_props["output_file_path"] = "prod_align_protein_ca.pdb"
    input_yaml_path = Path("pymol_align_protein_ca.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)
    assert Path("prod_align_protein_ca.pdb").exists()
