1. Open **SPGI**
2. Add a **Line chart**
3. Set **Split** = `Series` — the legend must list every distinct value of `Series` (typically 4–10 categories), each rendered with a distinct color from the categorical palette (no two adjacent categories share the same color)
4. On the **Context Panel > Data**, enable **Multi Axis** (`lc.props.multiAxis = true`) — each Y line gets its own subplot, and each subplot has its own legend driven by `Split`
5. Save layout → re-apply layout — `Split`, `Multi Axis`, and per-line legends must survive the round-trip
6. Save project → reopen project — same state must persist
7. Configure two Y columns: `lc.props.yColumnNames = ['Average Mass', 'TPSA']` — both lines render, each with its own legend block
8. Replace one Y column via the in-plot column selector (or `lc.props.yColumnNames = ['Average Mass', 'NIBR logP']` if the in-plot selector is hover-only) — the corresponding legend block must update to reflect the new column's values
9. Save layout → re-apply layout — the new Y column choice must persist
10. Save project → reopen project — same state must persist

---
{
  "order": 6,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
