#### Сhem filter

1. Open **spgi-100** 
2. Open the Filter Panel
3. In the **Structure** filter, draw the CCC(N(C)C)=O molecule. 
4. Hover over the **Structure** filter - the gear icon appears - click the gear icon - the searchTypeInput dropdown input appears. 
5. Sequentially switch between the dropdown values  and check the filter state:
  * Included in - 0 rows
  * Exact - 0 rows
  * Similar - 0 rows
  * Not contains - 86 rows
  * Not included in - 100 rows
  * Contains- 14 rows
6. Close All

#### Bio

1. Open peptides.csv
2. Wait for the AlignedSequence column to render sequences
3. Open the **Filter Panel**
4. For the **AlignedSequence** filter, in the Substructure input field, enter `T-T-Y-K-N-Y-V` - verify the result
5. Close All
---
{
"order": 2,
"datasets": ["System:AppData/Chem/tests/spgi-100.csv","System:DemoFiles/bio/peptides.csv"]
}
