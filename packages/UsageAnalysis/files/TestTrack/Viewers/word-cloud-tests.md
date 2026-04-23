# Word Cloud tests

Word Cloud is part of the Charts package (`public/packages/Charts`).

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Word Cloud (from Charts package)

## Column assignment

1. Open Properties, set Column to RACE — words appear for each race value, sized by frequency
2. Change Column to DIS_POP — word cloud updates to show disease population values
3. Change Column to SITE — word cloud updates to site values
4. Change Column to SEX — only two words appear (M/F)

## Shape

1. Open Properties, set Shape to circle — words arranged in a circular outline
2. Change Shape to diamond — words arranged in a diamond outline
3. Change Shape to triangle-forward — triangular arrangement
4. Change Shape to triangle
5. Change Shape to pentagon
6. Change Shape to star — star-shaped arrangement

## Font size range

1. Open Properties, set Min Text Size to 10 and Max Text Size to 60
2. Verify the most frequent words are larger and less frequent words are smaller
3. Set Min Text Size to 20 and Max Text Size to 20 — all words render at the same size
4. Restore Min Text Size to 12, Max Text Size to 48

## Text rotation

1. Open Properties, set Min Rotation Degree to 0 and Max Rotation Degree to 0 — all words are horizontal
2. Set Min Rotation Degree to -90 and Max Rotation Degree to 90 — words appear at various angles
3. Set Rotation Step to 45 — rotation angles are multiples of 45 degrees
4. Restore defaults

## Grid size

1. Open Properties, set Grid Size to 2 — words are packed more tightly
2. Set Grid Size to 20 — words are spaced further apart, fewer words fit
3. Restore Grid Size to 8

## Draw out of bound

1. Open Properties, check Draw Out Of Bound
2. Verify words can extend beyond the viewer border
3. Uncheck Draw Out Of Bound — words are clipped to the viewer bounds

## Font family

1. Open Properties, set Font Family to serif — words render in a serif font
2. Change Font Family to monospace — words render in monospace
3. Change Font Family back to sans-serif

## Bold

1. Open Properties, check Bold — all words become bold
2. Uncheck Bold — words return to normal weight

## Tooltip on hover

1. Hover over a word in the cloud — tooltip appears showing the word and its row count
2. Move to another word — tooltip updates

## Word click and selection

1. Click a word — corresponding rows are selected in the linked grid
2. Click another word — selection updates to new word's rows
3. Click empty space — selection is cleared

## Filter interaction

1. Open the filter panel, add a filter: SEX = M
2. Verify the word cloud updates — counts reflect only male patients, word sizes change
3. Remove the filter — word cloud restores original proportions

## Viewer resize

1. Resize the viewer by dragging its border — word cloud reflows to fill the new dimensions
2. Expand to full screen (Alt+F) — word cloud fills the screen
3. Press Alt+F again — viewer returns to normal size

## Context menu

1. Right-click on the word cloud — context menu appears
2. Verify standard viewer options are present (Properties, Save as PNG, etc.)
3. Click Properties — property panel opens

## Error state (too many categories)

1. Open Properties, set Column to a column with more than 500 unique values (e.g., a high-cardinality numerical column cast to string)
2. Verify the viewer shows an error message instead of the word cloud
3. Set Column back to RACE — error clears and word cloud renders normally

---
{
  "order": 100,
  "datasets": ["System:DemoFiles/demog.csv"]
}
