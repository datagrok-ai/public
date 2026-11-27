# Review: 03-more-viewers (More viewers)

Fresh-eyes review (2026-11-22, ~21:00-21:10): feature lines, descriptions.md, and a contact sheet at the key lines and
results of each clip, taken from the newest guide takes (the coordinator re-rendered them with the slower outlines,
so most are 10-25 s longer than the lengths in descriptions.md; the headings' m:ss need updating once the re-renders
are delivered). Judged on content, text and jumps, not pacing.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Compare groups with a box plot | OK. |
| 02 | Show shares with a pie chart | OK. |
| 03 | Small multiples with a trellis plot | OK. |
| 04 | Show a hierarchy with a tree map or a sunburst | OK, one jump: the tree map is docked to the right off camera (a mid-scenario `Given ... docks`), so it appears already moved. |
| 05 | Find related columns with a correlation plot | OK. |
| 06 | Compare many columns with a PC plot | OK. |
| 07 | Summarize every column with the Statistics viewer | OK, one jump: the viewer is docked to the bottom off camera. |
| 08 | Show flows with a Sankey diagram | OK. |
| 09 | Compare rows across many columns with a radar | OK. The second car's shape is drawn grey under the hover line (product bug 8 of the maker), so "Another car, another shape" is weaker than the first click. |
| 10 | Show rows as cards with the Tile Viewer and the Form | OK with two jumps: the Tile Viewer is docked off camera, and the context panel closes off camera before the Form is added. Also, on "Click a card: the form shows its row" the outline first lights row 1's card, then the lane scrolls and the press lands on row 3's card (the form does show row 3). |
| 11 | See where points crowd with a density plot | OK. |
| 12 | See a whole table at once with a heat map | OK. The grid beside it scrolls sideways (its model column goes out of view) once the heat map is sorted; harmless. |
| 13 | See every pair of columns with a matrix plot | OK. |
| 14 | Show links between things with a network diagram (1:01) | refilmed, still weak: the click on a character was written `And` after a `Then`, so it was filmed as a hidden check (no outline, 1 s); it is now a `When`, moved after the settings (a selection takes the context panel off the viewer's settings), and the clip ends on "33 selected rows" with the tooltip "Target: Eddard Stark". Still: the nodes carry no names (Node1 / Node2 Label set to Source / Target drew nothing, so that was dropped), and the physics layout keeps moving. A readable version needs labels that draw (product question) or a smaller table. |
| 15 | See dates on a calendar | OK. The calendar draws very large in the docked area (only Dec 1989 and Jan 1990 show); readable. |

Verdicts: 14 OK, 1 refilmed but still weak (14). Jumps to fix in a later pass: the toolbox is hidden off camera in the
11 clips that add their viewer from it, and 04, 07 and 10 dock a viewer or close the context panel off camera.

## Titles and descriptions

Plain and accurate. Every description starts "A <viewer> shows..." and then "This video adds one from...": consistent,
a little formulaic across 15 clips but not AI-sounding; left as is. No banned words.

## Systematic

- Jumps (the lead's 21:00 rule): every clip that adds a viewer from the toolbox then hides the toolbox with a
  mid-scenario `Given the toolbox pane is hidden`, so the toolbox vanishes between two steps; 04, 07 and 10 also dock
  the new viewer off camera. A filmed close of the toolbox (its tab click) or simply leaving it open would remove the
  jump; docking can be a filmed drag (`user docks ... viewer to ...` as a `When`).
- Sankey (08) and network diagram (14) open their table with Browse steps written as `Given`: fine as a starting
  state, since nothing was on screen before.

## No-jumps pass (2026-11-22, ~21:35-21:50)

Every feature now passes film.sh's no-jumps and hidden-gesture checks with no warning. Each changed clip was refilmed
in a page of its own, checked at the old jump with frames.sh, and delivered over the old video.

| NN | Clip | What changed | New length |
|---|---|---|---|
| 01 | Box plot | the toolbox is closed on camera with its sidebar icon, line "Hide the toolbox to give the plot room" | 1:22 |
| 02 | Pie chart | same, "...give the chart room" | 1:00 |
| 03 | Trellis plot | same, "...give the plot room" | 1:17 |
| 04 | Tree map / sunburst | toolbox closed on camera, then the map (which the toolbox docks under the grid) dragged by its title bar onto the dock wheel's right edge, line "Drag its title bar to dock it on the right"; claim `docked along the right edge` | 1:36 |
| 05 | Correlation plot | toolbox closed on camera | 1:06 |
| 06 | PC plot | toolbox closed on camera | 1:02 |
| 07 | Statistics | toolbox closed on camera, then the viewer dragged by its title bar to the bottom edge, line "Drag its title bar to the bottom, so every statistic fits"; claim `docked along the bottom edge` | 0:51 |
| 10 | Tile Viewer and Form | the viewer dragged by its title bar to the bottom on camera; the context panel hidden with the status bar's Context Panel toggle ("Hide the context panel for more room") before the Form is added. The card clicked for the Form is now row 2's, first in the None lane: row 3's sat half hidden under row 1's in the High lane, so the outline lit row 1's card while the lane scrolled under the press (the reviewer's note); now outline and press are on the same card, form reads 2, then 3 | 1:06 |
| 11 | Density plot | toolbox closed on camera | 0:55 |
| 12 | Heat map | toolbox closed on camera, "...give the map room" | 1:00 |
| 13 | Matrix plot | toolbox closed on camera; refilmed in its own page, the stray "↺" from 4.5 s is gone (checked at 4.6, 33, 45, 55, 70 s) | 1:12 |
| 15 | Calendar | the toolbox now stays open for the whole clip: closed (tried on camera), it let the calendar grow so large that January's label had no room and the month click lit a blank strip; open, the Jan and December labels are drawn and outlined | 0:37 |

- 09 (radar), checked for the "↺" at 55.5-61.5 s: refilmed once in its own page and it comes back. It is the product's
  own reset icon beside the Show Values row the pointer has just checked; it stays there for ~4 s after the pointer
  moves to the grid, until the next click repaints the panel. Not a cross-clip leftover, not floating: the delivered
  take is kept (the refilm is the same, 1:10).
- 08, 09 and 14 had no jump; their headings' m:ss now match the delivered files (0:43, 1:10, 1:00).
