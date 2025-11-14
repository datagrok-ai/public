#name: ReadSas
#language: python
#input: blob fileInput
#output: dataframe df

import pandas as pd
from io import BytesIO

# .xpt files are SAS Transport (XPORT) format
try:
    bytes_io_object = BytesIO(fileInput)
    df = pd.read_sas(bytes_io_object, format='xport', encoding='utf-8')
            
except Exception as e:
    raise Exception(f"Error reading SAS file: {str(e)}")