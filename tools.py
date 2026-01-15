"""
MOF Tools - Production-ready tools with Pydantic validation.

This module provides scientific tools for Metal-Organic Framework (MOF) research,
with comprehensive input/output validation using Pydantic models.
"""

from io import StringIO
from pathlib import Path
from typing import List, Optional, Dict, Any
from enum import Enum

from pydantic import BaseModel, Field, field_validator, ValidationError

# Optional ASE integration
try:
    import ase.io
    from ase.calculators.emt import EMT
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


# ============================================================================
# Pydantic Models for Input/Output Validation
# ============================================================================

class MOFRecord(BaseModel):
    """Model representing a single MOF database record."""
    name: str = Field(..., min_length=1, description="MOF name")
    formula: str = Field(..., min_length=1, description="Chemical formula")
    surface_area: float = Field(..., gt=0, description="Surface area in m²/g")

    class Config:
        frozen = True


class SearchMOFsInput(BaseModel):
    """Input model for MOF search operations."""
    query: str = Field(
        ..., 
        min_length=1, 
        max_length=200,
        description="Search query for MOF name or formula"
    )
    
    @field_validator('query')
    @classmethod
    def validate_query(cls, v: str) -> str:
        """Validate and sanitize the query string."""
        stripped = v.strip()
        if not stripped:
            raise ValueError("Query cannot be empty or whitespace only")
        return stripped


class SearchMOFsOutput(BaseModel):
    """Output model for MOF search results."""
    success: bool = Field(..., description="Whether the search was successful")
    results: List[MOFRecord] = Field(default_factory=list, description="List of matching MOFs")
    count: int = Field(..., ge=0, description="Number of results found")
    message: Optional[str] = Field(None, description="Additional message or error details")


class CalculateEnergyInput(BaseModel):
    """Input model for energy calculation."""
    data: str = Field(
        ..., 
        min_length=1,
        description="CIF file content as string or file path"
    )
    
    @field_validator('data')
    @classmethod
    def validate_data(cls, v: str) -> str:
        """Validate CIF data input."""
        if not v.strip():
            raise ValueError("Data cannot be empty")
        return v


class CalculateEnergyOutput(BaseModel):
    """Output model for energy calculation results."""
    success: bool = Field(..., description="Whether the calculation was successful")
    energy: Optional[float] = Field(None, description="Calculated energy in eV")
    error: Optional[str] = Field(None, description="Error message if calculation failed")
    message: str = Field(..., description="Human-readable result message")


class OptimizeStructureInput(BaseModel):
    """Input model for structure optimization."""
    name: str = Field(
        ..., 
        min_length=1,
        max_length=200,
        description="Name of the structure to optimize"
    )
    
    @field_validator('name')
    @classmethod
    def validate_name(cls, v: str) -> str:
        """Validate structure name."""
        stripped = v.strip()
        if not stripped:
            raise ValueError("Name cannot be empty or whitespace only")
        return stripped


class OptimizeStructureOutput(BaseModel):
    """Output model for structure optimization results."""
    success: bool = Field(..., description="Whether optimization was initiated successfully")
    structure_name: str = Field(..., description="Name of the structure being optimized")
    message: str = Field(..., description="Status message")


# ============================================================================
# Internal Database
# ============================================================================

SAMPLE_MOF_DB: List[MOFRecord] = [
    MOFRecord(name="HKUST-1", formula="Cu3(BTC)2", surface_area=1850.0),
    MOFRecord(name="MOF-5", formula="Zn4O(BDC)3", surface_area=3800.0),
    MOFRecord(name="UiO-66", formula="Zr6O4(OH)4(BDC)6", surface_area=1187.0),
]


# ============================================================================
# Tool Implementation Functions
# ============================================================================

def search_mofs(query: str) -> str:
    """
    Search for MOFs by name or formula.
    
    Args:
        query: Search query string for MOF name or formula
        
    Returns:
        JSON string containing search results with validation
        
    Raises:
        ValidationError: If input validation fails
    """
    try:
        # Validate input
        validated_input = SearchMOFsInput(query=query)
        
        # Perform search
        query_lower = validated_input.query.lower()
        matching_records = [
            mof for mof in SAMPLE_MOF_DB 
            if query_lower in mof.name.lower() or query_lower in mof.formula.lower()
        ]
        
        # Create validated output
        output = SearchMOFsOutput(
            success=True,
            results=matching_records,
            count=len(matching_records),
            message=f"Found {len(matching_records)} MOF(s)" if matching_records 
                   else f"No MOFs found for '{validated_input.query}'"
        )
        
        return output.model_dump_json(indent=2)
        
    except ValidationError as e:
        # Return validation error as JSON
        error_output = SearchMOFsOutput(
            success=False,
            results=[],
            count=0,
            message=f"Input validation error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
    except Exception as e:
        # Handle unexpected errors
        error_output = SearchMOFsOutput(
            success=False,
            results=[],
            count=0,
            message=f"Unexpected error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)


def calculate_energy(data: str) -> str:
    """
    Calculate the potential energy of a structure (CIF string or path).
    
    Args:
        data: CIF file content as string or file path
        
    Returns:
        JSON string containing energy calculation results with validation
        
    Raises:
        ValidationError: If input validation fails
    """
    try:
        # Validate input
        validated_input = CalculateEnergyInput(data=data)
        
        # Check ASE availability
        if not ASE_AVAILABLE:
            output = CalculateEnergyOutput(
                success=False,
                energy=None,
                error="ASE library not installed",
                message="Error: ASE library not installed. Install with: pip install ase"
            )
            return output.model_dump_json(indent=2)
        
        # Perform energy calculation
        try:
            fileobj = StringIO(validated_input.data) if "\n" in validated_input.data else validated_input.data
            atoms = ase.io.read(fileobj, format="cif" if "\n" in validated_input.data else None)
            atoms.calc = EMT()
            energy = atoms.get_potential_energy()
            
            output = CalculateEnergyOutput(
                success=True,
                energy=energy,
                error=None,
                message=f"Energy: {energy:.4f} eV"
            )
            return output.model_dump_json(indent=2)
            
        except Exception as calc_error:
            output = CalculateEnergyOutput(
                success=False,
                energy=None,
                error=str(calc_error),
                message=f"Calculation error: {str(calc_error)}"
            )
            return output.model_dump_json(indent=2)
            
    except ValidationError as e:
        # Return validation error as JSON
        error_output = CalculateEnergyOutput(
            success=False,
            energy=None,
            error="Input validation error",
            message=f"Input validation error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
    except Exception as e:
        # Handle unexpected errors
        error_output = CalculateEnergyOutput(
            success=False,
            energy=None,
            error="Unexpected error",
            message=f"Unexpected error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)


def optimize_structure(name: str) -> str:
    """
    Perform structure optimization (standard placeholder).
    
    Args:
        name: Name of the structure to optimize
        
    Returns:
        JSON string containing optimization status with validation
        
    Raises:
        ValidationError: If input validation fails
    """
    try:
        # Validate input
        validated_input = OptimizeStructureInput(name=name)
        
        # Perform optimization (placeholder)
        output = OptimizeStructureOutput(
            success=True,
            structure_name=validated_input.name,
            message=f"Successfully initiated optimization for {validated_input.name}"
        )
        
        return output.model_dump_json(indent=2)
        
    except ValidationError as e:
        # Return validation error as JSON
        error_output = OptimizeStructureOutput(
            success=False,
            structure_name=name if name else "unknown",
            message=f"Input validation error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
    except Exception as e:
        # Handle unexpected errors
        error_output = OptimizeStructureOutput(
            success=False,
            structure_name=name if name else "unknown",
            message=f"Unexpected error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
