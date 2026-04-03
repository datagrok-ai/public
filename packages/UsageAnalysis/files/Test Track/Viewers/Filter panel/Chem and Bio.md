#### Refresh with chem filter

1. Open **spgi-100** 
2. Open the Filter Panel
5. In the **Structure** filter, draw the CCC(N(C)C)=O molecule. 
6. Hover over the **Structure** filter - the gear icon appears - click the gear icon - the searchTypeInput dropdown input appears. 
6. Sequentially switch between the dropdown values  and check the filter state:
  * Included in - 0 rows
  * Exact - 0 rows
  * Similar - 0 rows
  * Not contains - 86 rows
  * Not included in - 100 rows
  * Contains- 14 rows
8. Go to the Toolbox, expand the File section and click **Refresh** — verify that filters are not cleared and data remains filtered (14 rows)
9. Close All

#### Bio

1. Open peptides.csv
1. Wait for the AlignedSequence column to render sequences
1. Open the **Filter Panel**
1. For the **AlignedSequence** filter, in the Substructure input field, enter `T-T-Y-K-N-Y-V` - verify the result
1. Close All
---
{
"order": 2,
"datasets": ["System:AppData/Chem/tests/spgi-100.csv","System:Demo/bio/peptides.csv"]
}
