---
feature: word-cloud
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [charts.viewer.word-cloud]
realized_as:
  - word-cloud-spec.ts
related_bugs: []
---

# Word Cloud (Playwright)

Word Cloud is part of the Charts package (`public/packages/Charts`).

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **Word Cloud** from **Toolbox > Viewers**

Everything below is driven through the UI: the Viewers toolbox, the Context Panel
property grid, and real mouse input on the canvas. The words that were actually
drawn (and their counts) are read back from the chart the viewer renders.

## Add the viewer

1. Click the **Word Cloud** icon in **Toolbox > Viewers**
2. Words are drawn, one per category of the selected column, and each word's
   count matches the number of rows with that value

## Column assignment

1. On **Context Panel > Data**, set **Column** to *RACE*
2. One word per race appears, each with the right row count, and the cloud is redrawn

## Too many categories

1. Set **Column** to *USUBJID* (5850 unique values)
2. The viewer shows "The Word cloud viewer requires categorical column with 500 or
   fewer unique categories" instead of a cloud
3. Set **Column** back to *RACE* — the error clears and the words come back

## Shape

1. On **Context Panel > Misc**, set **Shape** to *star*, then *diamond*, then *circle*
2. Each shape rearranges the words

## Font size

1. Set **Min Text Size** and **Max Text Size** both to *20* — the words are laid out
   at one size
2. Restore *12* and *48*

## Rotation

1. Set **Min Rotation Degree** and **Max Rotation Degree** to *0* — the words are
   laid out horizontally

## Bold and font family

1. Flip **Bold** — the words are restyled (it ships **on**)
2. Set **Font Family** to *monospace* — the words are restyled again
3. Restore both

## Tooltip

1. Hover a word — a tooltip appears showing that word's row count

## Selection

1. Click a word — exactly that word's rows are selected in the table

## Filtering

1. Filter the table to `SEX = M`
2. The viewer keeps drawing the same words without an error

   Note: the counts do **not** follow the filter — the cloud has no Row Source
   and keeps showing whole-table numbers.

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`word-cloud-spec.ts`. Kept verbatim so no scenario is lost.

### Font size range

> Manual

1. Open Properties, set Min Text Size to 10 and Max Text Size to 60
2. Verify the most frequent words are larger and less frequent words are smaller
3. Set Min Text Size to 20 and Max Text Size to 20 — all words render at the same size
4. Restore Min Text Size to 12, Max Text Size to 48

### Text rotation

> Manual

1. Open Properties, set Min Rotation Degree to 0 and Max Rotation Degree to 0 — all words are horizontal
2. Set Min Rotation Degree to -90 and Max Rotation Degree to 90 — words appear at various angles
3. Set Rotation Step to 45 — rotation angles are multiples of 45 degrees
4. Restore defaults

### Grid size

> Manual

1. Open Properties, set Grid Size to 2 — words are packed more tightly
2. Set Grid Size to 20 — words are spaced further apart, fewer words fit
3. Restore Grid Size to 8

### Draw out of bound

> Manual

1. Open Properties, check Draw Out Of Bound
2. Verify words can extend beyond the viewer border
3. Uncheck Draw Out Of Bound — words are clipped to the viewer bounds

### Viewer resize

> Manual

1. Resize the viewer by dragging its border — word cloud reflows to fill the new dimensions
2. Expand to full screen (Alt+F) — word cloud fills the screen
3. Press Alt+F again — viewer returns to normal size

### Context menu

> Manual

1. Right-click on the word cloud — context menu appears
2. Verify standard viewer options are present (Properties, Save as PNG, etc.)
3. Click Properties — property panel opens

---
{
  "order": 100,
  "datasets": ["System:DemoFiles/demog.csv"]
}
