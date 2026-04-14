1. Switching the bound table when multiple dataframes are open:
    1. Open SPGI, SPGI-linked1, SPGI-linked2.
    2. Open SPGI, add a **Form** viewer.
    3. Open viewer properties (**F4**) → **Data → Table**. Switch to SPGI-linked2, then SPGI-linked1, then back to SPGI. Verify fields rebind to the chosen table each time.
2. Building the form from a column list:
    1. On the Form viewer's ribbon, click the **list** icon — a column-picker dialog opens.
    2. Toggle several checkboxes; verify the form adds/removes the corresponding fields.
    3. Save a layout (**Ctrl+S** → Layout), reopen it, and verify the same fields are restored.
3. In design mode (ribbon **object-ungroup** icon on), drag a field's label to a new position and verify it stays there.
---
{
  "order": 9,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}