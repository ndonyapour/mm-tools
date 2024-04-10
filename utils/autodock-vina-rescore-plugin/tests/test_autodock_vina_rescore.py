"""Tests for autodock_vina_rescore."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_autodock_vina_rescore() -> None:
    """Test autodock_vina_rescore."""
    cwl_file = Path("autodock_vina_rescore.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)

    ligand_path = "ligand.pdbqt"
    ligand_path = str(Path(__file__).resolve().parent / Path(ligand_path))
    receptor_path = "receptor.pdbqt"
    receptor_path = str(Path(__file__).resolve().parent / Path(receptor_path))

    input_to_props["input_ligand_pdbqt_path"]["path"] = ligand_path
    input_to_props["input_receptor_pdbqt_path"]["path"] = receptor_path
    input_to_props["score_only"] = True
    input_to_props["output_log_path"] = "vina_rescore_pdbind.log"

    input_yaml_path = Path("autodock_vina_rescore.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)

    assert Path("vina_rescore_pdbind.log").exists()
