"""
Bandgap Predictor Tool

The bandgap predictor tool provides fast and scalable electronic bandgap prediction
for MOFs using a DPA-based property model, fine-tuned specifically for bandgap regression.
"""

from .base import (
    BaseModel, Field, ValidationError,
    Optional, Atoms, DeepProperty, os, DATA_DIR
)


class PredictBandgapInput(BaseModel):
    """Input model for bandgap prediction."""
    atoms_dict: dict = Field(
        ...,
        description="ASE Atoms object as dictionary (from parse_structure output)"
    )


class PredictBandgapOutput(BaseModel):
    """Output model for bandgap prediction results."""
    success: bool = Field(..., description="Whether the prediction was successful")
    bandgap: Optional[float] = Field(None, description="Predicted bandgap value in eV")
    error: Optional[str] = Field(None, description="Error message if prediction failed")
    message: str = Field(..., description="Human-readable result message")


def predict_bandgap(atoms_dict: dict) -> str:
    """
    Predict the electronic bandgap of a MOF structure using a DPA-based property model.
    
    Args:
        atoms_dict: ASE Atoms object as dictionary (from parse_structure output)
        
    Returns:
        JSON string containing bandgap prediction result
        
    Raises:
        ValidationError: If input validation fails
    """
    try:
        # Validate input
        validated_input = PredictBandgapInput(atoms_dict=atoms_dict)
        
        try:
            # Reconstruct Atoms object from dictionary
            atoms = Atoms(
                numbers=validated_input.atoms_dict["numbers"],
                positions=validated_input.atoms_dict["positions"],
                cell=validated_input.atoms_dict.get("cell"),
                pbc=validated_input.atoms_dict.get("pbc", [False, False, False])
            )
            #information for bandgap predictions
            coords = atoms.get_positions()
            cells = atoms.get_cell()
            atom_numbers = atoms.get_atomic_numbers()               
            atom_types = [x - 1 for x in atom_numbers]
            
            # The tool uses a DPA property model for bandgap regression
            model = DeepProperty(os.path.join(DATA_DIR, "bandgap.ckpt.pt"))
            
            # Predict bandgap value
            bandgap_val = model.eval(coords=coords, cells=cells, atom_types=atom_types)[0][0][0]
            
            output = PredictBandgapOutput(
                success=True,
                bandgap=float(bandgap_val),
                error=None,
                message=f"Bandgap prediction successful. Predicted value: {bandgap_val:.4f} eV"
            )
            return output.model_dump_json(indent=2)
            
        except Exception as calc_error:
            output = PredictBandgapOutput(
                success=False,
                bandgap=None,
                error=str(calc_error),
                message=f"Prediction error: {str(calc_error)}"
            )
            return output.model_dump_json(indent=2)
            
    except ValidationError as e:
        error_output = PredictBandgapOutput(
            success=False,
            bandgap=None,
            error="Input validation error",
            message=f"Input validation error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
    except Exception as e:
        error_output = PredictBandgapOutput(
            success=False,
            bandgap=None,
            error="Unexpected error",
            message=f"Unexpected error: {str(e)}"
        )
        return error_output.model_dump_json(indent=2)
