1. Open SPGI
1. Add a viewer. By default Row Source is set Filtered
1. Check Row Sources:
   1. Add some filters: on the Filter Panel and via the FIlter (${CAST Idea ID} <636500) setting on the Context Panel - viewer should respond to both filtrations
   1. Set Row source to All - viewer should respond to the FIlter setting  only
   1. Set Row source to Selected - viewer should be empty
   1. Reset filter on the Filter Panel and Select some rows - all selected rows should be displayed on the viewer (considering the Filter setting)
   1. Set Row source to SelectedOrCurrent - the same selected rows should be displayed on the viewer
   1. Reset selection and set any row as current - only current row should be displayed (considering the Filter setting)
   1. Set Row source to FIlteredSelected - the viewer should be empty
   1. Select some rows that intesect with filtered ones via the FIlter setting of the viewer  - only selected rows that correspond the Filter setting of the viewer should be displayed on the viewer. For example, data: SPGI,  FIlter =${CAST Idea ID} <636500 - select some rows where CAST Idea ID is from 636480 to 636504. Only rows from 636480 to 636499 should be displayed
   1. Set Row source to MouseOverGroup - the viewer is empty
   1. Add a pie chart (histogram, filter panel) and hover over its sections - the viewer should respond to hovering and the Filter setting
   1. Set Row source to CurrentRow - only current row that corresponds the Filter setting of the viewer should be displayed on the viewer.
   1. Set Row source to MouseOverRow - only row that corresponds the Filter setting of the viewer should be displayed on the viewer.
1. Open demog.csv
1. Set the Table setting of the viewer to demog, set Row Source to Filtered
1. Repeat step 3 using the demog (use Filter = ${AGE}>44)