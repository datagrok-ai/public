import csv
import random
import argparse
from typing import List, Dict, Callable, Any, Tuple

# --- Data Generation Functions ---
# Each function takes the plate number and the within-plate row number (0-95)
# and returns the corresponding value.

def get_well_layout_info(plate_number: int, row_number: int) -> Tuple[str, float, bool]:
    """
    Determines the well layout and returns (sample_id, concentration, is_control).
    
    Layout per row (12 wells):
    - Column 1: High Control (10.0 concentration)
    - Columns 2-9: Compound dilution series (8 concentrations)  
    - Columns 10-11: Empty/unused wells
    - Column 12: Low Control (0.0 concentration)
    
    Each plate row (A-H) tests a different compound.
    """
    plate_row = row_number // 12  # 0-7 (corresponds to rows A-H)
    plate_col = row_number % 12   # 0-11 (corresponds to columns 1-12)
    
    if plate_col == 0:  # Column 1 - High Control
        return "High Control", 10.0, True
    elif plate_col == 11:  # Column 12 - Low Control
        return "Low Control", 0.0, True
    elif 1 <= plate_col <= 8:  # Columns 2-9 - Compound dilution series
        # Calculate compound number: each row tests a different compound
        # Across multiple plates, continue numbering compounds
        compound_number = plate_row + 1 + (plate_number - 1) * 8
        sample_id = f"SMP-{compound_number}"
        
        # Standard 3-fold dilution series starting from 100
        dilution_step = plate_col - 1  # 0-7 for the 8 concentration points
        base_concentration = 100.0
        concentration = base_concentration / (3.0 ** dilution_step)
        # Round to reasonable precision
        concentration = round(concentration, 2)
        
        return sample_id, concentration, False
    else:  # Columns 10-11 - Empty wells
        return f"EMPTY-P{plate_number}-W{row_number + 1}", 0.0, False

def generate_barcode(plate_number: int, row_number: int) -> str:
    """Generates a unique barcode for each plate."""
    return f"LJ-PLATE-{plate_number}"

def generate_well_position(plate_number: int, row_number: int) -> str:
    """Calculates the well position (e.g., A1, H12) from a row index."""
    if not 0 <= row_number < 96:
        raise ValueError("Row number must be between 0 and 95.")

    row_letter = chr(ord('A') + (row_number // 12))
    column_number = (row_number % 12) + 1
    return f"{row_letter}{column_number}"

def generate_sample_id(plate_number: int, row_number: int) -> str:
    """Generates sample ID based on dose-response layout."""
    sample_id, _, _ = get_well_layout_info(plate_number, row_number)
    return sample_id

def generate_concentration(plate_number: int, row_number: int) -> float:
    """Generates concentration based on dose-response layout."""
    _, concentration, _ = get_well_layout_info(plate_number, row_number)
    return concentration

def generate_activity(plate_number: int, row_number: int) -> float:
    """
    Generates activity values that simulate dose-response behavior.
    For inhibition assays: higher concentration = lower activity
    Includes some random noise for realism.
    """
    sample_id, concentration, is_control = get_well_layout_info(plate_number, row_number)
    
    if is_control:
        if sample_id == "High Control":
            # High control should have low activity (strong inhibition)
            return round(random.uniform(0.05, 0.15), 3)
        else:  # Low Control
            # Low control should have high activity (no inhibition)
            return round(random.uniform(0.95, 1.05), 3)
    else:
        if "EMPTY" in sample_id:
            return 0.0
        
        # For compounds, simulate dose-response curve
        # Higher concentration = lower activity (inhibition)
        if concentration > 0:
            # Sigmoid-like response with noise
            # IC50 around concentration 3-10 for variability
            ic50 = random.uniform(3.0, 10.0)
            hill_slope = random.uniform(0.8, 1.5)
            
            # Hill equation: Activity = 1 / (1 + (concentration/IC50)^hill_slope)
            normalized_activity = 1.0 / (1.0 + (concentration / ic50) ** hill_slope)
            
            # Add some noise
            noise = random.uniform(-0.05, 0.05)
            activity = max(0.0, min(1.2, normalized_activity + noise))
            
            return round(activity, 3)
        else:
            return round(random.uniform(0.95, 1.05), 3)

def generate_random_float(plate_number: int, row_number: int) -> float:
    """Generates a random float between 0.0 and 1.0."""
    return round(random.random(), 3)

# --- Main Application Logic ---

def generate_plate_file(
    output_filename: str,
    column_names: List[str],
    num_plates: int
):
    """
    Generates a CSV file containing dose-response data for a specified number of 96-well plates.

    Args:
        output_filename: The name of the file to create (e.g., 'output.csv').
        column_names: A list of strings representing the desired column headers.
        num_plates: The number of plates to generate data for.
    """
    # A dictionary mapping known column names to their generator functions.
    generator_map: Dict[str, Callable[[int, int], Any]] = {
        'barcode': generate_barcode,
        'pos': generate_well_position,
        'SampleID': generate_sample_id,
        'Concentration': generate_concentration,
        'Activity': generate_activity,
    }

    rows_per_plate = 96  # Standard 8x12 plate

    print(f"Generating dose-response file '{output_filename}' with {num_plates} plate(s)...")
    print("Layout per plate:")
    print("  - Each row (A-H) tests a different compound")
    print("  - Column 1: High Control, Column 12: Low Control")
    print("  - Columns 2-9: 8-point dilution series per compound")
    print("  - Columns 10-11: Empty wells")
    print(f"  - Total compounds per plate: 8")
    print(f"  - Total compounds across all plates: {num_plates * 8}")

    try:
        with open(output_filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=column_names)
            writer.writeheader()

            # Loop through each plate
            for plate_idx in range(1, num_plates + 1):
                # Loop through each well in the plate
                for row_idx in range(rows_per_plate):
                    # Create a dictionary for the current row's data
                    row_data = {}
                    for col in column_names:
                        # Find the correct generator function for the column
                        generator_func = generator_map.get(col, generate_random_float)
                        row_data[col] = generator_func(plate_idx, row_idx)

                    writer.writerow(row_data)

        print("File generation complete.")
        print(f"Total rows written: {num_plates * rows_per_plate}")
        print(f"Dose-response series per compound: 8 concentrations")

    except IOError as e:
        print(f"Error: Could not write to file '{output_filename}'.")
        print(e)


def main():
    """
    Parses command-line arguments and runs the file generation logic.
    """
    parser = argparse.ArgumentParser(
        description="Generate mock CSV data for 96-well dose-response assay plates.\n"
                   "Each plate tests 8 compounds with 8-point dilution series plus controls.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '--columns',
        required=True,
        nargs='+',
        help="The list of column names for the CSV file.\n"
             "Example: --columns barcode pos SampleID Concentration Activity"
    )

    parser.add_argument(
        '--num_plates',
        required=True,
        type=int,
        help="The total number of plates to generate."
    )

    parser.add_argument(
        '--output',
        default='generated_drc_data.csv',
        type=str,
        help="The name of the output CSV file (default: generated_drc_data.csv)."
    )

    args = parser.parse_args()

    # Set random seed for reproducible results (optional)
    random.seed(42)

    generate_plate_file(args.output, args.columns, args.num_plates)


if __name__ == '__main__':
    main()
