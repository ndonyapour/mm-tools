"""Tests for gmx_rms_nofit."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_gmx_rms_nofit() -> None:
    """Test gmx_rms_nofit."""
    cwl_file = Path("gmx_rms_nofit.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)

    input_structure_path = "complex.gro"
    input_structure_path = str(
        Path(__file__).resolve().parent / Path(input_structure_path),
    )
    input_to_props["input_structure_path"]["path"] = input_structure_path

    input_traj_path = "prod_align_protein_ca.pdb"
    input_traj_path = str(Path(__file__).resolve().parent / Path(input_traj_path))
    input_to_props["input_traj_path"]["path"] = input_traj_path

    input_index_path = "MOL.ndx"
    input_index_path = str(Path(__file__).resolve().parent / Path(input_index_path))
    input_to_props["input_index_path"]["path"] = input_index_path

    input_to_props["output_xvg_path"] = "rmsd_equil_ligand_notime.xvg"
    input_yaml_path = Path("gmx_rms_nofit.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)
    assert Path("rmsd_equil_ligand_notime.xvg").exists()
