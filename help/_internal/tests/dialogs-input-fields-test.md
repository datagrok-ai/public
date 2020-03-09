<!-- TITLE: Tests: Dialogs Input Fields -->
<!-- SUBTITLE: -->

# Tests: Dialogs Input Fields

Most platform dialogs have input fields for data of different types. Testing the behavior of such places of the 
system when entering valid and non-valid ones is very important.

Testing the input fields of standard types (int, double, string, datetime) is reduced to checking the response of 
the system to input data of inappropriate types.

When inputting non-valid data into the input field, the following results are expected:
  
   * Field is displayed empty
  
   * Field will be highlighted in red
   
   * Tooltip about incorrect data type entered will appear
   
   * Further operation  will be stopped. Additional exceptions will not be thrown
   
When testing in this direction, it is necessary to check the behavior when inserting non-valid data types from the clipboard.

If you enter invalid data in dialog field, then after clicking ```Enter```, dialog should not be executed, and 
corresponding fields should be highlighted in red with a tooltip.

## Input Type "Column" 

Some dialog of platform use input fields with a column type.

For such fields there should be an additional sub-dialog ("Select Columns") under which you can select columns.
This sub-dialog contains a list of all available columns for selection, their data type, number of nulls, column search and column selector

When you open "Select Columns" sub-dialog, the following results are expected:

   * List shows all columns of the appropriate data type for the input field (all columns of the table, if there are no restrictions of the data type)
   
   * Displaying column data types
   
   * Number of nulls for each column is shown
   
   * Possibility to select certain columns, all columns and none
   
   * Columns search
   
   * When you change columns (deleting, renaming, changing the data type, changing the column data, etc.) changes are shown in the "Select columns" sub-dialog