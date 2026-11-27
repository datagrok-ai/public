# Review: 04-select-and-filter (Select and filter)

Fresh-eyes review (2026-11-22, ~21:05-21:20): feature lines, descriptions.md, and a contact sheet at the key lines and
results of each clip, from the takes the coordinator re-rendered with the slower outlines (delivered 21:05-21:07;
each video matches its guide take). The m:ss in the descriptions.md headings are the older, shorter renders and need
updating to the delivered lengths (e.g. 04 is now 1:03, 13 is 0:21). Judged on content, text and jumps.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Select rows | OK. |
| 02 | Filter rows with the filter panel | OK. |
| 03 | Categorical filter tricks | OK. |
| 04 | Filter with an expression | OK. |
| 05 | Filter with a hierarchy | OK. |
| 06 | Search text in a table | OK. Ctrl+F opens the toolbox's Search section on the left and squeezes the filter panel beside it; filmed as the key press, so not a jump. |
| 07 | Filter molecules by a substructure | OK. |
| 08 | Turn a selection into a filter | OK. |
| 09 | Filter from a viewer | OK. |
| 10 | Find and remove duplicate rows | OK. The status bar shows "Columns: 7" on iris after Selection to Top (the maker's product bug 3). |
| 11 | Filter to a list of IDs | text fixed: the description opened with a rhetorical question ("Got a list of compound IDs from a colleague?"); now a plain sentence. The pasted IDs show lowercased in the card's search box (product bug 2). |
| 12 | Save filters and apply them later | OK. |
| 13 | Copy the filtered rows into a new table | OK. The new view's top-right label still reads "72 filtered rows" (product bug 4). The dissolve into the new view doubles the grid for a moment around the Extract Rows press; transient. |
| 14 | Filter from a column header | OK. |
| 15 | Select rows with the same value | OK. |

Verdicts: 14 OK, 1 text fixed.

## Data claims (checked 2026-11-22 ~22:25, one frame of the delivered video per line that states a result)

All true, nothing changed. Every counted result is also asserted by the feature's own checks, and the frames show
the same numbers: 01: 7 selected in the status bar, 765 selected by the box on the plot and in the grid. 02: F keeps
3,243; turning the SEX card off goes from 2,078 to 3,558; reset gives all 5,850. 03: RA 2,550, UC 728, RA+UC 3,278,
inverted 2,572, "Ps" 1,409, a click on Psoriasis's count selects 1,209. 04: AGE > 50 keeps 2,176; AND gives 499;
the second rule off gives 2,176 again. 05: F keeps 3,243, the next level lists the four races, unticking Caucasian
leaves 420. 06: "chocolate" 25 (highlighted), AND 12, Origin = Belgium 14 (all Belgian). 07: 33 of 100 contain the
ring (highlighted in each), Not contains 67. 08: 751 selected, 751 filtered, No HEIGHT true in those 751. 09: the RA
bar keeps 2,550, empty space gives all, the zoom 3,962, Reset View all. 10: rows 102 and 143 are the same flower
(5.8, 2.7, 5.1, 1.9, virginica), 149 rows after removing one. 11: 5 compounds. 12: the preset brings back 37 Asian F
rows. 13: the new table has 72 rows, all Asian. 14: 1,889 rows after Esc. 15: 354 rows with RACE = Other; AGE > 70
selects 171 with nothing filtered.

## Titles and descriptions

Plain and accurate; Mac shortcuts are given (Cmd+click). One rhetorical question removed (11). The section avoids
off-camera jumps well: set-up `Given` steps only build the starting state (cards in the filter panel, a narrowed
grid), and every move after the first visible step is a filmed gesture.

## Systematic

- Several clips start with the filter panel already open and set up; that's a starting state, fine under the
  "no jumps" rule.
- The re-rendered dissolves sometimes catch a menu and the next view half-blended (13 at 14 s, 01 at 25 s); worth a
  look in the renderer if the lead finds it distracting.

## Merges (2026-11-22, ~23:15-23:35)

The lead found clips under ~35 s too thin. Two pairs merged into one clip each, filmed again, reviewed with the
checklist (pass, audit 0 off, timeline, frames at every note and result) and delivered under the kept number:

| NN | Clip | Took in | Length | Notes |
|---|---|---|---|---|
| 11 | Filter to a list of IDs and copy the rows into a new table | 13 (Extract Rows) | 0:46 | spgi: the + card for Id, the pasted list (5 rows), the Filtered count's menu, Extract Rows, then a click on the original's tab shows it still has 100 rows. First take: a `the video says` line right before a captioned menu pick showed ~1 s; made a caption on the click instead. The new view's top-right label reads "5 filtered rows" (product bug 4). |
| 14 | Select and filter from a column | 15 (Shift+Enter) | 0:45 | demog: Shift+Enter on RACE = Other (354), Esc, the Search box's Shift+Enter on AGE > 70 (171), then the AGE header's popup filter (1,889) and Add filter. The 171 selected rows show in orange in the popup's histogram; the last Esc (closing the popup) also clears the selection, which is what Esc does. |

Data claims: every count is asserted by the feature and visible in the status bar of the frame it is spoken over.
