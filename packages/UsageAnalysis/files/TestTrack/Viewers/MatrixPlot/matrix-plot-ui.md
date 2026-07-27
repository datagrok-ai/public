# Matrix plot manual checklist

Human-only visual checks for the Matrix plot; automated coverage lives in the
other scenarios in this folder. What remains here needs a human eye: the
Auto Layout resize heuristics that hide and restore labels and axes as the
viewer size changes.

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

---
{
  "order": 120,
  "datasets": ["System:DemoFiles/demog.csv"]
}
