"""
Base module for MOF tools.

Provides common imports, constants, and base classes used across all tools.
"""

from io import StringIO
from typing import List, Optional

from pydantic import BaseModel, Field, field_validator, ValidationError, ConfigDict

import ase.io
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS, LBFGS, FIRE
from ase.filters import FrechetCellFilter
try:
    from ase.constraints import FixSymmetry
except ImportError:
    # Some ASE versions might have it elsewhere or not at all
    FixSymmetry = None

__all__ = ['Atoms', 'EMT', 'BFGS', 'LBFGS', 'FIRE', 'FrechetCellFilter', 'FixSymmetry',
           'BaseModel', 'Field', 'field_validator', 'ValidationError', 
           'ConfigDict', 'Optional', 'List', 'StringIO', 'ase']
