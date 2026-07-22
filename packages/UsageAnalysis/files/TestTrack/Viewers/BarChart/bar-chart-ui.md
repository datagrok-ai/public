# Bar chart tests (manual)

Human-only visual/gesture checks. The automatable steps (basic click-to-select,
stacked-segment filtering, On Click modes, Show Selected Rows, double-click Reset View,
grid color-coding driving the bar colors with a layout round-trip, and context-menu
composition) have moved into the spec(s); what remains here needs a human eye —
modifier-key gestures, zoom/pan by mouse, hover tooltips, and eyeballing exact
color-scheme edits on the canvas.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Bar chart

## Additive selection gestures

1. Set **On Click** to **Select**
2. **Ctrl + click** a bar — selection adds that category
3. **Shift + drag** a rectangle over multiple bars — all covered categories are selected

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
