#name: PythonEnvNumpyCompat
#description: Verifies numpy ABI compatibility when pip packages depend on numpy
#language: python
#environment: env_numpy_compat
#meta.timeout: 600
#input: int rows
#output: int result

import numpy as np
import pandas as pd
import pyarrow as pa
import cv2

# The functional checks below are the real ABI test: a pip-installed opencv
# compiled against a different numpy major than the conda one fails cv2 import
# or the pyarrow interop with _ARRAY_API errors. (An old `numpy < 2` assert
# guarded the pyarrow-13 era base env; pyarrow >= 25 is numpy-2-native.)

# Verify pyarrow works (would fail with _ARRAY_API error if numpy ABI mismatch)
table = pa.table({'col': pa.array(range(rows))})

# Verify pandas+pyarrow interop
df = pd.DataFrame({'x': range(rows), 'y': range(rows)})
arrow_table = pa.Table.from_pandas(df)
df2 = arrow_table.to_pandas()

# Verify opencv works with the constrained numpy
img = np.zeros((10, 10, 3), dtype=np.uint8)
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

result = len(df2) + gray.shape[0]
