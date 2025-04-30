1. Open SPGI, SPGI-linked1, SPGI-linked2
1. Go to  Data > Link Tables
1. Set Tables to SPGI and SPGI-linked1
1. Set Link type to  'selection to Filter' 
1. Set  Key columns to  Id <> Concept Id
1. Click  Link
1. Click New Link
1. Set Tables to SPGI-linked2 and SPGI-linked1
1. Set Link type to  'filter to filter' 
1. Set Key columns to: 
   1. Sample Name
   1. link column 1
   1. link column 2
   1. link column 3
1.  Click Link and close the dialog
1. Go to **SPGI** and select some rows (from the beginning of the table)
1. Go to **SPGI-linked2**
1. Open the Filter Panel and filter out some rows
**Expected Result**: The SPGI-linked1 table should show the intersection of the rows selected in SPGI and filtered in SPGI-linked2.