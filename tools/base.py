"""
Base module for MOF tools.

Provides common imports, constants, and base classes used across all tools.
"""

from io import StringIO
from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, Field, field_validator, ValidationError, ConfigDict

import ase.io
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS, LBFGS, FIRE

__all__ = ['Atoms', 'EMT', 'BFGS', 'LBFGS', 'FIRE', 
           'BaseModel', 'Field', 'field_validator', 'ValidationError', 
           'ConfigDict', 'Optional', 'List', 'StringIO', 'ase']
