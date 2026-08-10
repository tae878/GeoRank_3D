"""GeoRank_3D: reliability-aware discontinuity extraction for rockfall source zones.

Reference implementation accompanying the audit of the Chilas slope manuscript.
See ``AUDIT_REPORT.md`` for the correspondence between each module and the
manuscript's equations and sections, and for the defects in the original
MATLAB release that this package corrects.
"""

__version__ = "1.0.0"

from .config import Config, ConfigError

__all__ = ["Config", "ConfigError", "__version__"]
