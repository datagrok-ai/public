1. Open the demog dataset.
2. Create a new calculated column (*Weight2*)
- Open the Add new column dialog.
- Enter the formula: *${WEIGHT} + 100*.
- Name the column *Weight2*.
- Press OK.
- Expected result: A new column Weight2 is created with values correctly calculated.
3. Create a new calculate column (*Weight3*)
- Open the Add new column dialog.
- Enter the formula: *${Weight2} + 100*.
- Name the column *Weight3*.
- Press OK.
- Expected result: A new column Weight3 is created with values correctly calculated based on Weight2.
4. Create a new calculate column (*Weight4*)
- Open the Add new column dialog.
- Enter the formula: *Log10(${Weight3}) - 0.2*.
- Name the column *Weight4*.
- Press OK.
- Expected result: A new column Weight4 is created with values correctly calculated based on Weight3.
5. Verify formula dependencies and recalculations
- Open the Context panel.
- Modify the formulas for Weight2, Weight3, and Weight4.
- Observe if the changes correctly propagate through all dependent columns.

Expected result:
- Changing the formula of Weight2 should automatically update Weight3 and Weight4.
- Changing the formula of Weight3 should update Weight4 accordingly.
- The recalculations should reflect the correct dependency hierarchy without errors.

Additional Notes:
- Ensure that modifying the formulas does not cause calculation errors or unexpected behavior.
- Validate that changes are saved and persist after closing and reopening the dataset.


---
{
  "order": 7
}