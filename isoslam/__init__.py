"""IsoSLAM."""

from importlib.metadata import version

release = version("isoslam")
__version__ = ".".join(release.split("."[:2]))
