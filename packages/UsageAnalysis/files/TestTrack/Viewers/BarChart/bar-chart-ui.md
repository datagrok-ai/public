# Bar chart tests (manual)

Automated coverage lives in the paired spec(s) of this section; this file lists
only checks that need a human eye.

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
> round-trip) is now covered automatically in barchart-setup-interact-spec.ts, and
> context-menu composition is now covered automatically in bar-chart-spec.ts.
> Eyeballing exact edited colors on the canvas stays a human-side check.

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv,System:AppData/Chem/tests/spgi-100.csv"]
}
