---
feature: barchart
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Automated coverage lives in the paired spec(s) of this section; this file
  lists only checks that need a human eye.
---

# Bar chart tests (manual)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Bar chart

## Zoom and pan

1. **Alt + drag** vertically on the chart — chart should zoom into that range of categories
2. Use **mouse wheel** to scroll through categories when zoomed in
3. **Drag** on empty space (no modifiers) — chart should pan vertically through categories

## Scrolling (range sliders)

1. Change the vertical range slider and use it to scroll
2. Set **Value** to CAST Idea ID (open spgi-100 dataset). Scroll through categories using the range slider

## Tooltip configuration

1. Right-click on the bar chart and open the **Tooltip** tab
2. Verify that all tooltip options are functioning correctly
3. Hover over different bars — tooltip content should match configuration

> Grid color coding driving the bar colors (with a save → close → reload layout
> round-trip) and context-menu composition are covered automatically.
> Eyeballing exact edited colors on the canvas stays a human-side check.

## Pick Up / Apply

1. Customize the bar chart: set Category to RACE, Stack to SEX, and change the color scheme
2. Add a second Bar chart with default settings
3. Right-click the customized bar chart and select Pick Up / Apply > Pick Up
4. Right-click the second bar chart and select Pick Up / Apply > Apply — the second viewer now matches the first (category, stacking, colors)

## Cleanup

1. Close all — the scenario only adds in-session viewers; nothing persists on the server.

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
