"""
MOF Tools Package

Production-ready tools for Metal-Organic Framework (MOF) research with
comprehensive input/output validation using Pydantic models.

This package provides three main tools:
- parse_structure: Load and validate structure formats (CIF, XYZ, POSCAR)
- static_calculation: Perform static energy evaluation using DPA models
- optimize_geometry: Perform structure relaxation using DPA force fields
"""

from .parse_structure import parse_structure
from .static_calculation import static_calculation
from .optimize_geometry import optimize_geometry
from .predict_bandgap import predict_bandgap
from .fetch_structure import fetch_structure

__all__ = [
    'parse_structure',
    'static_calculation',
    'optimize_geometry',
    'predict_bandgap',
    'fetch_structure',
]

__version__ = '1.0.0'
