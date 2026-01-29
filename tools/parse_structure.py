"""
Structure Parser Tool

Load and validate different structure formats (CIF, XYZ, POSCAR) into ASE Atoms objects.
"""

from .base import (
    BaseModel, Field, field_validator, ValidationError,
    Optional, StringIO, ase, ASE_AVAILABLE
)


class ParseStructureInput(BaseModel):
    """Input model for structure parsing."""
    data: str = Field(
        ..., 
        min_length=1,
        description="Structure file content as string or file path"
    )
    format: Optional[str] = Field(
        None,
        description="Structure format (cif, xyz, vasp). Auto-detected if not provided"
    )
    
    @field_validator('data')
    @classmethod
    def validate_data(cls, v: str) -> str:
        """Validate structure data input."""
        if not v.strip():
            raise ValueError("Data cannot be empty")
        return v
    
    @field_validator('format')
    @classmethod
    def validate_format(cls, v: Optional[str]) -> Optional[str]:
        """Validate format string."""
        if v is not None:
            allowed_formats = ['cif', 'xyz', 'vasp', 'poscar']
            v_lower = v.lower()
            if v_lower not in allowed_formats:
                raise ValueError(f"Format must be one of {allowed_formats}")
            return v_lower
        return v


class ParseStructureOutput(BaseModel):
    """Output model for structure parsing results."""
    success: bool = Field(..., description="Whether parsing was successful")
    atoms_dict: Optional[dict] = Field(None, description="ASE Atoms object as dictionary")
    num_atoms: Optional[int] = Field(None, description="Number of atoms in the structure")
    formula: Optional[str] = Field(None, description="Chemical formula")
    error: Optional[str] = Field(None, description="Error message if parsing failed")
    message: str = Field(..., description="Human-readable result message")


def parse_structure(data: str, format: Optional[str] = None) -> str:
    """
    Load and validate different structure formats into ASE Atoms object.
    
    Args:
        data: Structure file content as string or file path
        format: Structure format (cif, xyz, vasp). Auto-detected if not provided
        
    Returns:
        JSON string containing parsed structure with validation
        
    Raises:
        ValidationError: If input validation fails
    """
    try:
        # Validate input
        validated_input = ParseStructureInput(data=data, format=format)
        
        # Check ASE availability
        if not ASE_AVAILABLE:
            output = ParseStructureOutput(
                success=False,
                atoms_dict=None,
                num_atoms=None,
                formula=None,
                error="ASE library not installed",
                message="Error: ASE library not installed. Install with: pip install ase"
            )
            return output.model_dump_json(indent=2)
        
        # Parse structure
        try:
            # Determine if data is file content or path
            is_file_content = "\n" in validated_input.data or len(validated_input.data) > 500
            
            if is_file_content:
                fileobj = StringIO(validated_input.data)
                file_format = validated_input.format if validated_input.format else "cif"
                atoms = ase.io.read(fileobj, format=file_format)
            else:
                # Assume it's a file path
                atoms = ase.io.read(validated_input.data, format=validated_input.format)
            
            # Convert Atoms object to dictionary for JSON serialization
            atoms_dict = {
                "positions": atoms.get_positions().tolist(),
                "numbers": atoms.get_atomic_numbers().tolist(),
                "cell": atoms.get_cell().tolist() if atoms.cell is not None else None,
                "pbc": atoms.get_pbc().tolist() if atoms.pbc is not None else [False, False, False],
            }
            
            output = ParseStructureOutput(
                success=True,
                atoms_dict=atoms_dict,
                num_atoms=len(atoms),
                formula=atoms.get_chemical_formula(),
                error=None,
                message=f"Successfully parsed structure: {atoms.get_chemical_formula()} ({len(atoms)} atoms)"
            )
            return output.model_dump_json(indent=2)
            
        except Exception as parse_error:
            output = ParseStructureOutput(
                success=False,
                atoms_dict=None,
                num_atoms=None,
                formula=None,
                error=str(parse_error),
                message=f"Parsing error: {str(parse_error)}"
            )
            return output.model_dump_json(indent=2)
            
    except ValidationError as e:
        error_output = ParseStructureOutput(
            success=False,
            atoms_dict=None,
            num_atoms=None,
            formula=None,
            error="Input validation error",
            message=f"Input validation error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
    except Exception as e:
        error_output = ParseStructureOutput(
            success=False,
            atoms_dict=None,
            num_atoms=None,
            formula=None,
            error="Unexpected error",
            message=f"Unexpected error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
