#name: ReadSas
#language: python
#input: file fileInput
#output: dataframe df

import pandas as pd
from io import BytesIO

# .xpt files are SAS Transport (XPORT) format
try:
    df = pd.read_sas(fileInput, format='xport', encoding='utf-8')
            
except Exception as e:
    raise Exception(f"Error reading SAS file: {str(e)}")