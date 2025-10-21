### Add sticky metadata to a single cell

1. Open SPGI dataset.
2. Select a row, click on cell under column Structure. 
3. Open Context Panel → Sticky Meta. Find TestSchema1 (created in previous test). Fill metadata:
   - set rating = 5, 
   - notes = "test note", 
   - verified = true
   - review_date = current date
4. Save the schema.
- Expected Result:
  - Blue marker (circle) appears in top-right of the cell (with a tooltip on its hover);
  - Metadata values are visible in Context Panel.
  - After closing and reopening the dataset, same metadata is still present.

### Create sticky column

1. In SPGI table, via Context Panel or column menu, add sticky column to the analysis.
2. Verify column appears with header marker (circle).
3. Check that rows with metadata show proper values.
4. Sort/filter by sticky column.
5. Remove sticky column from view (without deleting metadata).
- Expected Result:
   - Sticky column behaves like a regular column: sorting/filtering works; removing the column does not delete metadata.
   - After removal and re-adding sticky column, metadata values reappear.

### Batch edit metadata on multiple rows

1. Select multiple rows in SPGI (e.g. via shift-click)
2. In popup menu choose Sticky Meta > Edot for all properties (Rows = Selected)
3. Set verified = true, notes = "batch note"
- Expected Result:
   - Metadata applied to all selected rows.
   - All selected rows show the new metadata in tooltips; no missing rows.

---
{
"order": 2,
"datasets": ["System:DemoFiles/SPGI.csv"]
}