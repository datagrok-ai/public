#name: Revvity ETL
#description: ETL process for Revvity Signals data - exports libraries, processes CSV files, and prepares data for MolTrack
#language: python
#output: dataframe result [Processed dataframe with all libraries combined]

import pandas as pd
import requests
import time
import re
from typing import Dict, List, Optional, Any, Tuple

# API configuration
API_URL = 'https://dgrok-snb-intg.signalsresearch2.revvitycloud.com/api/rest/v1.0'
API_KEY = 'GoBhcedLSUHiXwzu0a772JwAB1Qp5hAVglAQ+rqpfmMuMBPDi92kQLKWFyv+6O5J2bZELg=='

HEADERS = {
    'x-api-key': API_KEY,
    'Accept': 'application/vnd.api+json',
}

# Type mappings
NON_STRING_TYPES_MAPPING = {
    'DECIMAL': 'double',
    'PERCENTAGE': 'double',
    'MASS': 'double',
    'MOLECULAR_MASS': 'double',
    'DATETIME': 'datetime',
    'VOLUME': 'double',
    'COUNT': 'int',
}

TYPES_TO_PARSE_FOR_UNITS = ['PERCENTAGE', 'MASS', 'MOLECULAR_MASS', 'VOLUME']


def get_substring_before_last_dash(s: str) -> str:
    """Extract substring before the last dash."""
    last_dash_index = s.rfind('-')
    if last_dash_index == -1:
        return s
    return s[:last_dash_index]


def find_first_non_null_val(series: pd.Series) -> Optional[str]:
    """Find first non-null value in a pandas Series."""
    for val in series:
        if pd.notna(val):
            return str(val)
    return None


def extract_value_and_units(s: str) -> Optional[Dict[str, str]]:
    """Extract numeric value and units from a string."""
    if not s or pd.isna(s):
        return None
    s = str(s)
    # Regex to match number (with commas) and units
    regex = r'([+-]?(?:\d{1,3}(?:,\d{3})+|\d*\.?\d+))\s*([a-zA-Z%°\'"²³µ/]*)$'
    match = re.match(regex, s)
    
    if match:
        return {
            'value': match.group(1).replace(',', ''),  # Remove commas for numeric conversion
            'units': match.group(2) if match.group(2) else ''
        }
    return None


def update_prop_with_units(props: List[Dict], col_name: str, units: str):
    """Update property units in the properties list."""
    for prop in props:
        if prop['name'] == col_name:
            prop['units'] = units
            break


def create_property(field: Dict, calculated: bool = False) -> Dict:
    """Create a property definition from a field."""
    prop = {
        'name': field['name'],
        'value_type': NON_STRING_TYPES_MAPPING.get(field.get('dataType', ''), 'string')
    }
    if field.get('calculated') == True or calculated:
        prop['property_class'] = 'calculated'
    return prop


def get_libraries() -> List[Dict]:
    """Fetch libraries from the API."""
    # Note: Original TypeScript code has double slash, but single slash should work
    url = f'{API_URL}/materials/libraries'
    response = requests.get(url, headers=HEADERS)
    response.raise_for_status()
    libraries = response.json()
    libraries_data = libraries.get('data', [])
    if not isinstance(libraries_data, list):
        libraries_data = [libraries_data] if libraries_data else []
    return libraries_data


def get_file(lib_name: str, start_after_id: Optional[str] = None) -> Tuple[pd.DataFrame, int]:
    """
    Request bulk export file and wait for completion, then download and return as DataFrame.
    Returns (dataframe, total_count).
    """
    # Request bulk export
    url = f'{API_URL}/materials/{lib_name}/bulkExport?fileType=CSV&chemFormat=SMILES&limit=25000'
    if start_after_id:
        url += f'&startAfter={start_after_id}'
    
    response = requests.post(url, headers=HEADERS)
    response.raise_for_status()
    file_info = response.json()
    file_id = file_info['data']['id']
    report_id = file_info['data']['attributes']['reportId']
    
    # Poll for completion
    status = ''
    total_count = 0
    
    print(f'Waiting for bulk export to complete (report ID: {report_id})...')
    while status != 'COMPLETED':
        time.sleep(5)  # Wait 5 seconds between polls
        report_url = f'{API_URL}/materials/bulkExport/reports/{report_id}'
        report_response = requests.get(report_url, headers=HEADERS)
        report_response.raise_for_status()
        report_status = report_response.json()
        
        if not total_count:
            total_count = report_status['data']['attributes'].get('total', 0)
        
        status = report_status['data']['attributes']['status']
        count = report_status['data']['attributes'].get('count', 0)
        
        if total_count > 0:
            progress = (count / total_count) * 100
            print(f'Progress: {progress:.1f}% - Loaded {count} of {total_count}')
        else:
            print(f'Status: {status}, Count: {count}')
    
    print('Export completed. Downloading file...')
    
    # Download the file
    download_url = f'{API_URL}/materials/bulkExport/download/{file_id}'
    download_headers = {
        'x-api-key': API_KEY,
        'Accept': '*/*',
    }
    file_response = requests.get(download_url, headers=download_headers)
    file_response.raise_for_status()
    
    # Parse CSV
    csv_content = file_response.text
    df = pd.read_csv(pd.StringIO(csv_content))
    
    return df, total_count


def process_library(lib: Dict) -> Tuple[pd.DataFrame, List[Dict], List[Dict]]:
    """Process a single library: export data, transform columns, and return properties."""
    lib_name = lib['attributes']['name']
    print(f'\nProcessing library: {lib_name}')
    
    # Create property definitions
    columns_to_parse = []
    
    compounds_props = []
    if 'assets' in lib['attributes'] and 'fields' in lib['attributes']['assets']:
        compounds_props = [create_property(it) for it in lib['attributes']['assets']['fields']]
        for prop in compounds_props:
            if prop['value_type'] in TYPES_TO_PARSE_FOR_UNITS:
                columns_to_parse.append(prop['name'])
    
    if 'assetCalculatedFields' in lib['attributes']:
        calculated_props = [create_property(it, True) for it in lib['attributes']['assetCalculatedFields']]
        compounds_props.extend(calculated_props)
        for prop in calculated_props:
            if prop['value_type'] in TYPES_TO_PARSE_FOR_UNITS:
                columns_to_parse.append(prop['name'])
    
    batch_props = []
    if 'batches' in lib['attributes'] and 'fields' in lib['attributes']['batches']:
        batch_props = [create_property(it) for it in lib['attributes']['batches']['fields']]
        for prop in batch_props:
            if prop['value_type'] in TYPES_TO_PARSE_FOR_UNITS:
                columns_to_parse.append(prop['name'])
    
    if 'batchCalculatedFields' in lib['attributes']:
        calculated_batch_props = [create_property(it, True) for it in lib['attributes']['batchCalculatedFields']]
        batch_props.extend(calculated_batch_props)
        for prop in calculated_batch_props:
            if prop['value_type'] in TYPES_TO_PARSE_FOR_UNITS:
                columns_to_parse.append(prop['name'])
    
    # Add revvity corporate id properties
    compounds_props.append({'name': 'revvityCorporateCompoundId', 'value_type': 'string'})
    batch_props.append({'name': 'revvityCorporateBatchId', 'value_type': 'string'})
    
    # Get first file
    total_df, total_count = get_file(lib_name)
    print(f'Initial file loaded: {len(total_df)} rows, total expected: {total_count}')
    
    # Get additional files if needed (pagination)
    while len(total_df) < total_count:
        last_id = total_df['ID'].iloc[-1]
        print(f'Fetching next batch starting after ID: {last_id}')
        res_df, _ = get_file(lib_name, last_id)
        total_df = pd.concat([total_df, res_df], ignore_index=True)
        print(f'Total rows so far: {len(total_df)} of {total_count}')
    
    # Rename ID column to revvityCorporateBatchId and create revvityCorporateCompoundId column
    if 'ID' in total_df.columns:
        total_df = total_df.rename(columns={'ID': 'revvityCorporateBatchId'})
        total_df['revvityCorporateCompoundId'] = total_df['revvityCorporateBatchId'].apply(get_substring_before_last_dash)
    
    # Rename structure column
    if 'Chemical Structure (SMILES)' in total_df.columns:
        total_df = total_df.rename(columns={'Chemical Structure (SMILES)': 'smiles'})
    
    # Parse units and extract numeric values for columns that need it
    for col_name in columns_to_parse:
        if col_name in total_df.columns:
            col = total_df[col_name]
            # Find first non-null value to extract units
            first_val = find_first_non_null_val(col)
            if first_val:
                res = extract_value_and_units(first_val)
                if res:
                    update_prop_with_units(compounds_props, col_name, res['units'])
                    update_prop_with_units(batch_props, col_name, res['units'])
            
            # Create new column with extracted numeric values
            new_col_name = f'{col_name}_2'
            def extract_numeric_value(x):
                if pd.isna(x):
                    return None
                result = extract_value_and_units(str(x))
                if result:
                    try:
                        return float(result['value'])
                    except (ValueError, TypeError):
                        return None
                return None
            
            total_df[new_col_name] = col.apply(extract_numeric_value)
            # Remove old column and rename new one
            total_df = total_df.drop(columns=[col_name])
            total_df = total_df.rename(columns={new_col_name: col_name})
    
    print(f'Library {lib_name} processed: {len(total_df)} rows')
    print(f'Compounds properties: {len(compounds_props)}')
    print(f'Batch properties: {len(batch_props)}')
    
    # TODO: update MolTrack schema with compoundsProps, batchProps; load totalDf to MolTrack
    
    return total_df, compounds_props, batch_props


# Main execution
print('Starting Revvity ETL process...')

# Get libraries
libraries = get_libraries()
print(f'Found {len(libraries)} libraries')

# Process each library
for lib in libraries:
    try:
        df, compounds_props, batch_props = process_library(lib)
        print(f'\nCompounds properties: {compounds_props}')
        print(f'Batch properties: {batch_props}')
        print(f'Dataframe: {len(df)} rows')
    except Exception as e:
        print(f'Error processing library {lib.get("attributes", {}).get("name", "unknown")}: {e}')
        import traceback
        traceback.print_exc()

# Return empty dataframe as placeholder
result = pd.DataFrame()
result

