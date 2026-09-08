"""
Identifier:     KSC-SJA-__init__.py
Name:           __init__.py
Description:    FGS guider matching tool of CSST
Author:         Huimei Feng
Created:        2024-01-05
Modified-History:
    2024-01-05, Huimei Feng, created
"""

from .Ecliptic_polyMatch import core_poly_match

__version__ = "0.0.1"
__all__ = ["core_poly_match"]
