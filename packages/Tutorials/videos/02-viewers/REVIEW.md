# Review: 02-viewers (Viewers)

Fresh-eyes review of the delivered clips (2026-11-22, ~21:00): feature lines, descriptions.md, and a contact sheet at
the key lines and results of each clip. Clips 02-05 were checked from the then-delivered .mp4s. After this review the
coordinator refilmed 03, 04 and 05 (their 20:38-20:40 takes carried the stray "↺" property-reset icon) and re-rendered
every clip with the slower outlines; the content judged here is unchanged, the lengths in descriptions.md are not.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Add viewers to a table (0:54) | OK. Both ways in (toolbox, gallery with a tag and a search), docking, full screen and close are all on screen. |
| 02 | Configure a scatter plot (1:10) | OK. Zoom shows "11 filtered rows", the regression line, square markers, labels, X histogram and r = 0.826 all land. |
| 03 | Configure a histogram (0:50) | OK. |
| 04 | Configure a bar chart (0:55) | OK. |
| 05 | Configure a line chart (0:51) | OK. The med pick shows "med" on the axis; the product bug where the list later blanks is not in this clip. |
| 06 | Viewers work together (0:32) | OK. |
| 07 | Change a viewer's settings (0:46) | OK. |
| 08 | Summarize a table with a pivot table (0:37) | OK. |
| 09 | Choose what a tooltip shows (0:30) | OK. While unticking, each column's own tooltip pops up over the plot; a little busy, not wrong. |
| 10 | Copy a viewer or its look (0:26) | OK. |
| 11 | Show only some rows in a viewer (0:39) | caption fixed in the data-claim pass (re-rendered, no refilm): "The table and the other viewers keep all rows" showed on a view with no other viewer (only the grid); now "The grid keeps all 1,000 rows". The plot stays empty for about 8 s after Row Source = Selected, before rows are selected; it reads as intended. |
| 12 | Use the legend (0:32) | OK. No stray icon (the round icon beside a legend item is the palette, shown on hover). |

Verdicts: 12 OK.

## Data claims (checked 2026-11-22 ~22:20, one frame of the delivered video per line that states a result)

All true; one line reworded. 02: the zoom leaves "11 filtered rows", the fit shows r = 0.826. 03: the hovered bin
holds 792 rows (the same 792 Show Values writes over it), the range handles filter to 5319, the double-click brings
back all 5,850, the split draws F and M. 04: 6 DIS_POP bars, 12 stack segments, "by name" gives AS, Indigestion,
PsA, Psoriasis, RA, UC, then upright, then 100% bars. 05: months on X, 2 lines after the split, 2 charts, then 1
chart with 2 Y axes, med, spline. 06: a 597-row bin lights up on the scatter plot and the bar chart; 1,741 rows
selected by the Shift+drag; 2,550 by the RA bar; the Psoriasis filter leaves 1,209 everywhere. 08: 24 rows by
DIS_POP and RACE, 4 by RACE, F/M columns after the pivot, the new table has 4 rows. 09: the tooltip lists exactly the
7 ticked columns ("7 of 11 checked"). 10: the second plot takes the hollow markers and both histograms and keeps AGE
against HEIGHT; Clone makes a third. 11: 466 rows over 44 on the plot; the 8 of the 25 selected rows that are over
44 (counted from the grid's ages). 12: one category shows 59 rows, the table keeps 1,000, Ctrl+click adds Black.

Fixed: 11's line "The table and the other viewers keep all rows" showed on a view with no other viewer; it is now
"The grid keeps all 1,000 rows", which the status bar shows (re-rendered from the same take, delivered at 0:39). The
description now says "the grid and any other viewer keep all rows". 01 and 07 state nothing about data.

## Titles and descriptions

Plain and accurate; no banned words; lengths in the headings match the videos. No text changes made.

## Systematic

- No jumps (the lead's 21:00 rule): no feature of this section has a `Given` after its first visible step; the toolbox
  is closed with a filmed click on its tab (01).

- The section reads well as a set: each viewer clip goes columns on the plot, then a right-click, then the gear and
  one or two settings per section, as the lead asked.
- "Ctrl+click" is shown with a Ctrl key cap; the descriptions add "Cmd+click on a Mac" where it matters (12). Keep that
  pattern in other sections.
