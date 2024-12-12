"""Test the utilities module and functions."""

from pathlib import Path

import pytest

from isoslam import utils


@pytest.mark.parametrize(
    ("sample_config", "new_values", "expected_config"),
    [
        pytest.param(
            {
                "a": 1,
                "b": 2,
                "c": "something",
                "base_dir": "here",
                "output_dir": "there",
            },
            {"c": "something new"},
            {
                "a": 1,
                "b": 2,
                "c": "something new",
                "base_dir": Path("here"),
                "output_dir": Path("there"),
            },
            id="Change value of c",
        ),
        pytest.param(
            {
                "a": 1,
                "b": 2,
                "c": "something",
                "base_dir": "here",
                "output_dir": "there",
                "nested_dictionary": {"e": 1},
            },
            {"e": 1000000},
            {
                "a": 1,
                "b": 2,
                "c": "something",
                "base_dir": Path("here"),
                "output_dir": Path("there"),
                "nested_dictionary": {"e": 1000000},
            },
            id="Change nested value e",
        ),
    ],
)
def test_update_config(sample_config: dict, new_values: dict, expected_config: dict) -> None:
    """Test updating configuration."""
    updated_config = utils.update_config(sample_config, new_values)

    assert isinstance(updated_config, dict)
    assert updated_config == expected_config
    assert isinstance(updated_config["base_dir"], Path)
    assert isinstance(updated_config["output_dir"], Path)
