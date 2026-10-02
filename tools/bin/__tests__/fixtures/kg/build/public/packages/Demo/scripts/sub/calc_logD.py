#name: Calculate logD
#description: Calculates logD.
#language: python
#environment: demo-env
#reference: https://en.wikipedia.org/wiki/Partition_coefficient
#sample: chem/smiles.csv
#tags: demo
#input: dataframe table
#input: column molecules {semType: Molecule}
#input: double pH = 7.4 {caption: pH}
#output: dataframe result {action:join(table)}
#test: expect(1, 1) //cat: Types
#meta.timeout: 900000

import pandas as pd
