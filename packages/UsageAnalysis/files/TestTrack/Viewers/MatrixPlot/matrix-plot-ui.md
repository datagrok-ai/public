---
feature: matrixplot
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  What remains here needs a human eye: the Auto Layout resize heuristics that
  hide and restore labels and axes as the viewer size changes, and
  cross-viewer point interactions inside cells.
---

# Matrix plot manual checklist

Human-only visual checks for the Matrix plot

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Matrix plot

## Auto Layout — Resize Behavior

1. Open the viewer settings. Verify **Auto Layout** (Style section) is true by default.
2. Resize the viewer to a small size.
3. Verify column labels and axes hide automatically.
4. Resize the viewer to a large size.
5. Verify column labels and axes reappear.

## Interact with Plot Elements

(Close all, open spgi-100, and add a Matrix plot for this scenario)

1. Hover over a data point in an off-diagonal cell — the corresponding row is
   highlighted in the grid
2. Click a data point — the corresponding row becomes current in the grid
3. Drag across a region in a cell — all data points in the region are selected
   in the grid

---
{
  "order": 120,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
