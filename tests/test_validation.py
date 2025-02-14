"""Test validation function."""

from collections.abc import Callable
from contextlib import nullcontext as does_not_raise
from pathlib import Path

import pytest
from schema import Or, Schema, SchemaError

from isoslam import validation

TEST_SCHEMA = Schema(
    {
        "a": Path,
        "b": Or("aa", "bb", error="Invalid value in config, valid values are 'aa' or 'bb"),
        "logical": bool,
        "positive_integer": lambda n: 0 < n,
        "int_or_float": Or(int, float, error="Invalid value in config should be type int or float"),
    }
)


@pytest.mark.parametrize(
    ("config", "expectation"),
    [
        pytest.param(
            {"a": Path(), "b": "aa", "logical": True, "positive_integer": 4, "int_or_float": 10.0},
            does_not_raise(),
            id="valid configuration",
        ),
        pytest.param(
            {"a": "path", "b": "aa", "logical": True, "positive_integer": 4, "int_or_float": 10.0},
            pytest.raises(SchemaError),
            id="invalid Path",
        ),
        pytest.param(
            {"a": Path(), "b": 3, "logical": False, "positive_integer": 4, "int_or_float": 10.0},
            pytest.raises(SchemaError),
            id="invalid str",
        ),
        pytest.param(
            {"a": Path(), "b": "bb", "logical": True, "positive_integer": -4, "int_or_float": 10.0},
            pytest.raises(SchemaError),
            id="negative integer",
        ),
        pytest.param(
            {"a": Path(), "b": "aa", "logical": False, "positive_integer": 4, "int_or_float": "five"},
            pytest.raises(SchemaError),
            id="invalid int/float",
        ),
        pytest.param(
            {"a": Path(), "b": "bb", "logical": "True", "positive_integer": 4, "int_or_float": "five"},
            pytest.raises(SchemaError),
            id="invalid boolean",
        ),
        pytest.param(
            {"a": Path(), "b": "bb", "logical": 1, "positive_integer": 4, "int_or_float": 2},
            pytest.raises(SchemaError),
            id="boolean as int (1)",
        ),
        pytest.param(
            {"a": Path(), "b": "aa", "logical": 0, "positive_integer": 4, "int_or_float": 4},
            pytest.raises(SchemaError),
            id="boolean as int (0)",
        ),
    ],
)
def test_validate_config(config: dict, expectation: Callable) -> None:
    """Test various configurations."""
    with expectation:
        validation.validate_config(config, schema=TEST_SCHEMA, config_type="Test YAML")
