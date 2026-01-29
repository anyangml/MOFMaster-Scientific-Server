"""
Base module for MOF tools.

Provides common imports, constants, and base classes used across all tools.
"""

from io import StringIO
from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, Field, field_validator, ValidationError, ConfigDict

# Optional ASE integration
try:
    import ase.io
    from ase import Atoms
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS, LBFGS, FIRE
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    Atoms = None

__all__ = ['ASE_AVAILABLE', 'Atoms', 'EMT', 'BFGS', 'LBFGS', 'FIRE', 
           'BaseModel', 'Field', 'field_validator', 'ValidationError', 
           'ConfigDict', 'Optional', 'List', 'StringIO', 'ase']
