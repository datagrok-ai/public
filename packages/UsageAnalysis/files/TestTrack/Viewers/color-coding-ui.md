### Prerequisites
demog table is open with no color coding applied.

### Grid-level color coding via Context Menu
1. Right-click the grid header → **Grid Color Coding → All**.
   Expected: every column shows color coding.

2. Right-click the grid header → **Grid Color Coding → Color scheme for table** → pick a different scheme.
   Expected: all colored columns update to the selected scheme.

3. Right-click the grid header → **Grid Color Coding → None**.
   Expected: all color coding is removed.

4. Right-click the grid header → **Grid Color Coding → Auto**.
   Expected: color coding is restored to the state before step 3.

### Layout save and restore
*(Precondition: demog is open with the custom color coding from color-coding.md step 2.)*

5. Add any viewer (e.g., Scatter Plot) to the table view.

6. Save the view layout (toolbar or context menu).

7. Close the demog table. Reopen it from **Files**, apply the saved layout, then close the added viewer.
   Expected: the table shows the same color coding state as step 2.

### Linked color coding — save and reload
*(Precondition: demog open with linked color coding from step 5.)*

- Save layout; reload → verify RACE, HEIGHT, and chain columns remain Linked with correct sources.
- Save project; reopen → verify the same.

### Edit Color Scheme dialog (SPGI_v2 dataset)
8. Load **SPGI_v2**; apply **Linear** color coding to the `CAST Idea ID` column.

9. Right-click the column header → **Color Coding → Edit...** → click **Edit scheme...**.

10. In the dialog:
    - Hover over **"Use absolute values"** → verify a tooltip appears.
    - Check the checkbox → verify absolute value input fields become active.
    - Modify min, max, and inflection point colors → verify the live preview in the grid updates immediately.
    - Click **Cancel** → verify all in-dialog changes are reverted.
    - Click **+** to add a color range entry → verify the scheme preview updates.
    - Click **−** to remove the added entry → verify the scheme reverts.

11. Click **OK** → verify the final scheme is applied to the grid.

12. Click the **arrows icon** (invert) next to the scheme dropdown → verify the gradient is reversed in the grid.
