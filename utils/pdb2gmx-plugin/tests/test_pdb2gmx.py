"""Tests for pdb2gmx."""
import sys
from pathlib import Path

current_dir = Path(__file__).resolve().parent
target_dir = current_dir.parent.parent.parent / "cwl_utils"
sys.path.append(str(target_dir))

from cwl_utilities import call_cwltool  # noqa: E402
from cwl_utilities import create_input_yaml  # noqa: E402
from cwl_utilities import parse_cwl_arguments  # noqa: E402


def test_pdb2gmx() -> None:
    """Test pdb2gmx."""
    cwl_file = Path("pdb2gmx.cwl")
    input_to_props = parse_cwl_arguments(cwl_file)
    input_pdb_path = "protein_1e3g.pdb"
    input_pdb_path = str(Path(__file__).resolve().parent / Path(input_pdb_path))
    input_to_props["input_pdb_path"]["path"] = input_pdb_path
    input_to_props[
        "config"
    ] = '{"water_type":"spce","force_field":"amber99sb-ildn","ignh":true,"merge":false}'
    input_to_props["output_gro_path "] = "receptor.gro"
    input_to_props["output_top_zip_path"] = "system.zip"

    input_yaml_path = Path("pdb2gmx.yml")
    create_input_yaml(input_to_props, input_yaml_path)
    call_cwltool(cwl_file, input_yaml_path)

    assert Path("system.zip").exists()
