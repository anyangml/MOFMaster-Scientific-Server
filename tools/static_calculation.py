"""
Static Calculation Tool

Perform static energy evaluation using DPA models without modifying geometry.
Computes total energy, forces, and virial in a non-iterative manner.
"""

from .base import (
    BaseModel, Field, ValidationError,
    Optional, Atoms, DP, os, DATA_DIR
)


class StaticCalculationInput(BaseModel):
    """Input model for static energy calculation."""
    atoms_dict: dict = Field(
        ...,
        description="ASE Atoms object as dictionary (from parse_structure output)"
    )
    normalize_per_atom: bool = Field(
        False,
        description="Whether to compute energy per atom"
    )
    compute_forces: bool = Field(
        True,
        description="Whether to compute atomic forces"
    )
    compute_virial: bool = Field(
        False,
        description="Whether to compute virial tensor"
    )


class StaticCalculationOutput(BaseModel):
    """Output model for static energy calculation results."""
    success: bool = Field(..., description="Whether the calculation was successful")
    total_energy: Optional[float] = Field(None, description="Total energy in eV")
    energy_per_atom: Optional[float] = Field(None, description="Energy per atom in eV/atom")
    forces: Optional[list] = Field(None, description="Atomic forces in eV/Å, shape (N, 3)")
    virial: Optional[list] = Field(None, description="Virial tensor in eV, shape (3, 3)")
    error: Optional[str] = Field(None, description="Error message if calculation failed")
    message: str = Field(..., description="Human-readable result message")


def static_calculation(
    atoms_dict: dict,
    normalize_per_atom: bool = False,
    compute_forces: bool = True,
    compute_virial: bool = False
) -> dict:
    """
    Perform static energy evaluation using DPA model without modifying geometry.
    
    This is a non-iterative calculation that computes total energy, forces,
    and virial directly from the DPA model. The structure is not modified.
    
    Args:
        atoms_dict: ASE Atoms object as dictionary (from parse_structure output)
        normalize_per_atom: Whether to compute energy per atom
        compute_forces: Whether to compute atomic forces
        compute_virial: Whether to compute virial tensor
        
    Returns:
        Dictionary containing static calculation results with validation
        
    Raises:
        ValidationError: If input validation fails
    """
    try:
        # Validate input
        validated_input = StaticCalculationInput(
            atoms_dict=atoms_dict,
            normalize_per_atom=normalize_per_atom,
            compute_forces=compute_forces,
            compute_virial=compute_virial
        )
        
        # Perform static calculation
        try:
            # Reconstruct Atoms object from dictionary
            import numpy as np
            
            atoms = Atoms(
                numbers=validated_input.atoms_dict["numbers"],
                positions=validated_input.atoms_dict["positions"],
                cell=validated_input.atoms_dict.get("cell"),
                pbc=validated_input.atoms_dict.get("pbc", [False, False, False])
            )
            
            atoms.calc = DP(os.path.join(DATA_DIR, "forcefield.ckpt.pt"))
            
            # Compute total energy (static, non-iterative)
            total_energy = atoms.get_potential_energy()
            
            # Compute energy per atom if requested
            energy_per_atom = None
            if validated_input.normalize_per_atom:
                energy_per_atom = total_energy / len(atoms)
            
            
            # Compute forces if requested
            forces = None
            if validated_input.compute_forces:
                forces_array = atoms.get_forces()
                forces = forces_array.tolist()
            
            # Compute virial if requested
            virial = None
            if validated_input.compute_virial:
                try:
                    # Get stress tensor and convert to virial
                    stress = atoms.get_stress(voigt=False)  # 3x3 tensor
                    volume = atoms.get_volume()
                    virial_array = -stress * volume  # Virial = -stress * volume
                    virial = virial_array.tolist()
                except Exception:
                    # Some calculators don't support stress/virial
                    virial = None
            
            # Build message
            msg_parts = [f"Total energy: {total_energy:.4f} eV"]
            if energy_per_atom is not None:
                msg_parts.append(f"Energy/atom: {energy_per_atom:.4f} eV/atom")
            if forces is not None:
                max_force = np.max(np.abs(forces))
                msg_parts.append(f"Max force: {max_force:.4f} eV/Å")
            
            output = StaticCalculationOutput(
                success=True,
                total_energy=float(total_energy),
                energy_per_atom=float(energy_per_atom) if energy_per_atom is not None else None,
                forces=forces,
                virial=virial,
                error=None,
                message="Static calculation completed. " + ", ".join(msg_parts)
            )
            return output.model_dump()
            
        except Exception as calc_error:
            output = StaticCalculationOutput(
                success=False,
                total_energy=None,
                energy_per_atom=None,
                forces=None,
                virial=None,
                error=str(calc_error),
                message=f"Calculation error: {str(calc_error)}"
            )
            return output.model_dump()
            
    except ValidationError as e:
        # Return validation error as JSON
        error_output = StaticCalculationOutput(
            success=False,
            total_energy=None,
            energy_per_atom=None,
            forces=None,
            virial=None,
            error="Input validation error",
            message=f"Input validation error: {str(e)}"
        )
        return error_output.model_dump()
    except Exception as e:
        # Handle unexpected errors
        error_output = StaticCalculationOutput(
            success=False,
            total_energy=None,
            energy_per_atom=None,
            forces=None,
            virial=None,
            error="Unexpected error",
            message=f"Unexpected error: {str(e)}"
        )
        return error_output.model_dump()
