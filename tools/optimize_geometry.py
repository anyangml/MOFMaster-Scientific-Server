"""
Geometry Optimization Tool

Perform structure relaxation for MOFs using DPA machine-learning force fields.
"""

from .base import (
    BaseModel, Field, field_validator, ValidationError,
    Optional, Atoms, DP, BFGS, LBFGS, FIRE, FrechetCellFilter, FixSymmetry, os, DATA_DIR
)


class OptimizeGeometryInput(BaseModel):
    """Input model for geometry optimization."""
    atoms_dict: dict = Field(
        ...,
        description="ASE Atoms object as dictionary (from parse_structure output)"
    )
    fmax: float = Field(
        0.05,
        gt=0,
        description="Force convergence threshold in eV/Å"
    )
    max_steps: int = Field(
        200,
        gt=0,
        le=1000,
        description="Maximum number of optimization steps"
    )
    optimizer: str = Field(
        "BFGS",
        description="Optimizer type: BFGS, LBFGS, or FIRE"
    )
    relax_cell: bool = Field(
        True,
        description="Whether to relax lattice parameters"
    )
    fix_symmetry: bool = Field(
        True,
        description="Whether to maintain symmetry during optimization"
    )
    
    @field_validator('optimizer')
    @classmethod
    def validate_optimizer(cls, v: str) -> str:
        """Validate optimizer type."""
        allowed = ['BFGS', 'LBFGS', 'FIRE']
        v_upper = v.upper()
        if v_upper not in allowed:
            raise ValueError(f"Optimizer must be one of {allowed}")
        return v_upper


class OptimizationMetadata(BaseModel):
    """Metadata for optimization results."""
    converged: bool = Field(..., description="Whether optimization converged")
    final_fmax: float = Field(..., description="Final maximum force in eV/Å")
    steps: int = Field(..., description="Number of optimization steps performed")
    initial_energy: Optional[float] = Field(None, description="Initial energy in eV")
    final_energy: Optional[float] = Field(None, description="Final energy in eV")


class OptimizeGeometryOutput(BaseModel):
    """Output model for geometry optimization results."""
    success: bool = Field(..., description="Whether optimization completed successfully")
    optimized_atoms_dict: Optional[dict] = Field(None, description="Optimized ASE Atoms object as dictionary")
    metadata: Optional[OptimizationMetadata] = Field(None, description="Optimization metadata")
    error: Optional[str] = Field(None, description="Error message if optimization failed")
    message: str = Field(..., description="Human-readable result message")


def optimize_geometry(
    atoms_dict: dict,
    fmax: float = 0.05,
    max_steps: int = 200,
    optimizer: str = "BFGS",
    relax_cell: bool = False,
    fix_symmetry: bool = True
) -> dict:
    """
    Perform geometry optimization using DPA machine-learning force fields.
    
    Args:
        atoms_dict: ASE Atoms object as dictionary (from parse_structure output)
        fmax: Force convergence threshold in eV/Å
        max_steps: Maximum number of optimization steps
        optimizer: Optimizer type (BFGS, LBFGS, or FIRE)
        relax_cell: Whether to relax lattice parameters
        
    Returns:
        Dictionary containing optimization results with metadata
        
    Raises:
        ValidationError: If input validation fails
    """
    try:
        # Validate input
        validated_input = OptimizeGeometryInput(
            atoms_dict=atoms_dict,
            fmax=fmax,
            max_steps=max_steps,
            optimizer=optimizer,
            relax_cell=relax_cell,
            fix_symmetry=fix_symmetry
        )
        
        # Perform optimization
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
            
            # Get initial energy and forces
            initial_energy = atoms.get_potential_energy()
            
            # Select optimizer
            optimizer_class = {
                "BFGS": BFGS,
                "LBFGS": LBFGS,
                "FIRE": FIRE
            }[validated_input.optimizer]
            
            # Apply constraints and filters
            if validated_input.fix_symmetry and FixSymmetry is not None:
                atoms.set_constraint(FixSymmetry(atoms))
            
            # Use FrechetCellFilter for cell relaxation if requested
            target_atoms = atoms
            if validated_input.relax_cell:
                target_atoms = FrechetCellFilter(atoms)
            
            # Run optimization
            opt = optimizer_class(target_atoms)
            opt.run(fmax=validated_input.fmax, steps=validated_input.max_steps)
            
            # Get final results
            final_energy = atoms.get_potential_energy()
            final_forces = atoms.get_forces()
            final_fmax = np.max(np.abs(final_forces))
            converged = final_fmax <= validated_input.fmax
            
            # Convert optimized structure to dictionary
            optimized_dict = {
                "positions": atoms.get_positions().tolist(),
                "numbers": atoms.get_atomic_numbers().tolist(),
                "cell": atoms.get_cell().tolist() if atoms.cell is not None else None,
                "pbc": atoms.get_pbc().tolist() if atoms.pbc is not None else [False, False, False],
            }
            
            # Create metadata
            metadata = OptimizationMetadata(
                converged=converged,
                final_fmax=float(final_fmax),
                steps=opt.get_number_of_steps(),
                initial_energy=float(initial_energy),
                final_energy=float(final_energy)
            )
            
            output = OptimizeGeometryOutput(
                success=True,
                optimized_atoms_dict=optimized_dict,
                metadata=metadata,
                error=None,
                message=f"Optimization {'converged' if converged else 'did not converge'} after {opt.get_number_of_steps()} steps. "
                       f"Final fmax: {final_fmax:.4f} eV/Å, Energy: {final_energy:.4f} eV"
            )
            return output.model_dump()
            
        except Exception as opt_error:
            output = OptimizeGeometryOutput(
                success=False,
                optimized_atoms_dict=None,
                metadata=None,
                error=str(opt_error),
                message=f"Optimization error: {str(opt_error)}"
            )
            return output.model_dump()
            
    except ValidationError as e:
        error_output = OptimizeGeometryOutput(
            success=False,
            optimized_atoms_dict=None,
            metadata=None,
            error="Input validation error",
            message=f"Input validation error: {str(e)}"
        )
        return error_output.model_dump()
    except Exception as e:
        error_output = OptimizeGeometryOutput(
            success=False,
            optimized_atoms_dict=None,
            metadata=None,
            error="Unexpected error",
            message=f"Unexpected error: {str(e)}"
        )
        return error_output.model_dump()
