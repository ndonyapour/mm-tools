"""Tests for flatten_2d_to_1d."""
from pathlib import Path

from sophios.api.pythonapi import Step
from sophios.api.pythonapi import Workflow


def test_flatten_2d_to_1d() -> None:
    """Test flatten_2d_to_1d tool."""
    cwl_file_str = "flatten_2d_to_1d.cwl"
    cwl_file = Path(__file__).resolve().parent.parent / Path(cwl_file_str)
    flatten_2d_to_1d = Step(clt_path=cwl_file)
    flatten_2d_to_1d.input_2d_array = [["1"], ["2", "3"]]

    steps = [flatten_2d_to_1d]
    filename = "flatten_2d_to_1d"
    viz = Workflow(steps, filename)

    viz.run()
