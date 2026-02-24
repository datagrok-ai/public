#name: RowsTest
#language: python
#output: dataframe res
#test: CvmTests:getColumn(RowsTest(), 'calories').length == 3
import pandas as pd
data = {
  "calories": [420, 380, 390]
}
res = pd.DataFrame(data)
