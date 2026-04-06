### Linked Color-Coding Between Columns
Preconditions:
- Open the demog dataset.

#### Link RACE to WEIGHT (background coloring)
1. WEIGHT column:
   - Popup menu → Color Coding → switch to Categorical.
2. RACE column:
   - Popup menu → Color Coding → Linked.
   - Click "Edit…" and set Source column = WEIGHT.
Expected:
- RACE column background is colored using WEIGHT’s colors.

#### Link HEIGHT to WEIGHT (text coloring)
3. Select HEIGHT column.
4. In the Context Panel → Colors:
   - Type: Linked
   - Apply to: Text
   - Source column: WEIGHT
Expected:
- HEIGHT column text is colored using WEIGHT’s colors.

#### Source color changes propagate
5. Modify WEIGHT column color-coding:
   - Change colors
   - Switch between Linear and Conditional types
Expected:
- RACE and HEIGHT columns automatically update and inherit WEIGHT’s color changes.

####  Additional Checks
6. Save and reopen:
   - Save layout and reload it.
   - Save project and reopen it.
Expected:
- Linked color-coding works correctly after reload.

7. Linking chain (up to 5 levels):
   - Example: Column A → B → C → D → E
   - Each column linked to the previous one.
Expected:
- Up to 5 linked columns behave correctly and inherit colors through the chain.
