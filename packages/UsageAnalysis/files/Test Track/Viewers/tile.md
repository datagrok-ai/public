#### Tile viewer

1. **Multiple DataFrames Handling**
   - Open **SPGI**, **SPGI-linked1**, and **SPGI-linked2** datasets.
   - On the **SPGI** view, add a **Tile viewer**.
   - In the Context Panel, open **Data > Table** and switch between SPGI-linked2, SPGI-linked1, SPGI.
   - **Expected**: viewer rebinds to the chosen table and re-renders without errors.

2. **Building from List of Columns**
   - Right-click the viewer → **Edit form...**.
   - Drag a few columns onto the card, click **Close and apply**.
   - Save the layout (Ctrl+S) and reopen the project.
   - **Expected**: the saved card layout is restored exactly.

3. **Edit form with Table Change**
   - Open **demog** in addition to SPGI.
   - On SPGI, open **Edit form...**, change the source table to **demog**, press **Reset**, add a few demog columns, click **Close and apply**.
   - In the Tile viewer settings, switch **Data > Table** to **demog**.
   - **Expected**: viewer displays demog data; no errors and no red "broken-binding" placeholders inside tiles.

4. **Hamburger menu**
   - Open the Hamburger menu and walk through every item — each opens its dialog/submenu without errors.
   - Slowly hover **Properties > Data > Lanes**.
   - **Expected**: submenu opens smoothly; the page does not freeze even when the lanes column has many categories.

5. **Calculated column**
   - On SPGI, create a calculated column `${Average Mass} + 5` named `cc`.
   - Hamburger → **Edit form...** → **Reset**, add `cc` to the card, **Close and apply**.
   - Click the `cc` column header; in the Context Pane, change the formula and **Apply**.
   - **Expected**: column values and tiles update immediately.

6. **Calculated column with filters**
   - Apply filters on `cc`, **Stereo category**, and one more column.
   - Change the `cc` formula again.
   - **Expected**: tiles update correctly while filters stay applied.


---
{
  "order": 10
}
