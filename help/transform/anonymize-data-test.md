<!-- TITLE: Tests: Anonymize Data -->
<!-- SUBTITLE: -->

# Tests: Anonymize Data

Sometimes, it is needed to prepare a dataset that conveys an idea of the structure in the data along
with the patterns in it, but does not contain the real data points. "[Anonymize Data](anonymize-data.md)" functionality
can be handy here.


## Testing scenarios

1. Open "demog" dataset
  
1. Open "Anonymize Data..." dialog from  **Tools** menu

1. Select a numerical column "Age"

1. Click "OK" button
   * Values in the "Age" column changed with the randomization factor equal to 0.1 (
default value)

1. Anonymize data of string column "Sex"
   * Values ​​in the "Sex" column have changed to "Sex 0" and "Sex 1"
  
1. Use different values in "Number randomization factor" field (negatives, integers, floats) for anonymize data of "Height" field
   * Values ​​in the "Height" column are changed according to the randomization factor
   
1. Test non-functional modules (UI, popup menu, help, navigation, properties, etc.)
   * Non-functional modules work correctly and are intuitive

See also:
  * [Anonymize data](anonymize-data.md)
  * [Anonymize data Auto Test](anonymize-data-test.side)
