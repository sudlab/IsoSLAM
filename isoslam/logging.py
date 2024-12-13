"""Configure logging."""

import sys

from loguru import logger


def setup(level: str = "INFO") -> None:
    """
    Loguru setup with the required logging level and format.

    Parameters
    ----------
    level : str
        Log level, default is "INFO", other options "WARNING", "DEBUG" etc.
    """
    logger.remove()
    logger.add(sys.stderr)
    logger.add(
        sys.stderr,
        colorize=True,
        level=level.upper(),
        format="<green>{time:HH:mm:ss}</green> "
        "| <level>{level}</level> | "
        "<magenta>{file}</magenta>:<magenta>{module}</magenta>:<magenta>{function}</magenta>:<magenta>{line}</magenta>"
        " | <level>{message}</level>",
    )
