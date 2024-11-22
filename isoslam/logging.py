"""Configure logging."""

import sys

from loguru import logger

logger.remove()
logger.add(sys.stderr)
logger.add(
    sys.stderr,
    colorize=True,
    format="{time:HH:mm:ss} | <level>{level}</level> |<magenta>{file}</magenta>:<magenta>{module}</magenta>:"
    "<magenta>{function}</magenta>:<magenta>{line}</magenta> | <level>{message}</level>",
)
