---
feature: viewers
realizes_atlas: []   # untyped: cross-viewer docking, no single-feature atlas
priority: p2
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Dock targets are canvas overlays without stable drop targets for a script,
  so drag-docking is not automatable.
---

# Viewers docking — manual checklist

All scenarios start with:
1. Close all
2. Open demog
3. Add a Scatter plot, a Histogram, and a Bar chart

## Docking to view edges

1. Drag the Histogram by its header — docking overlay appears with drop zones
2. Drop it on the right edge zone of the view — Histogram docks to the right,
   full height; the other viewers shrink to make room
3. Drag the Bar chart to the bottom edge zone — it docks along the bottom,
   full width
4. Drag the splitter between two docked viewers — both resize accordingly,
   no rendering artifacts

## Docking relative to another viewer

1. Drag the Histogram over the Scatter plot — the overlay shows zones relative
   to that viewer
2. Drop on the Scatter plot's bottom zone — Histogram docks directly under the
   Scatter plot, sharing its width
3. Drop the Bar chart on the Scatter plot's left zone — it docks to the left
   of the Scatter plot only, not the whole view

## Layout persistence

1. Arrange the three viewers into a non-default docking (e.g. Scatter plot
   left, Histogram right-top, Bar chart right-bottom)
2. Save the layout (View > Layout > Save to Gallery)
3. Select View > Layout > Download — a `.layout` file is downloaded
4. Close All, open demog again
5. Apply the saved layout — all three viewers reappear in the same docking
   arrangement with the same relative sizes
6. Close All, open demog, drag the downloaded `.layout` file onto the view —
   the same layout is applied
7. Cleanup: delete the saved layout from the gallery and delete the
   downloaded file

---
{
  "order": 17,
  "datasets": ["System:DemoFiles/demog.csv"]
}
