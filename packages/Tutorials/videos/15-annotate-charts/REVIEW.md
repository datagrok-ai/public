# Review: 15-annotate-charts (Annotate and explain charts)

Fresh-eyes review (2026-11-22, ~23:00-23:10) of the 7 delivered clips: feature, descriptions.md block, and contact
sheets at the key lines and results of each delivered take (every delivered video matched its latest take's length
before the review). Judged on: is the one idea clear to someone new to Datagrok, is each line true when it shows,
plain wording, jumps. The m:ss in descriptions.md match the delivered videos.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Draw reference lines with formula lines | OK. Add Line opens the dialog at the median (77.4), typed to 100; the BMI 25 curve, its title and red color; Tools > Formula Lines... lists both. "Points above the curve have a BMI over 25" is true by the formula (544 of the 872 plotted rows, recomputed from data/demo/demog-1000.csv). |
| 02 | Mark a group of points with an annotation region | refilmed (0:44). The click selected the 26 rows, but in this plot the selection color is the same orange as versicolor: the only sign was the box's green points turning orange, with 2 s left in the clip. A closing hover on the status bar's "Selected: 26" now carries "The status bar counts them: 26 selected". The drag draws a visible box; the tooltip reads "Overlap / Versicolor and virginica mix here / 26 rows". |
| 03 | Label points with structures | OK. Structure and Id labels drawn while the column picker is open, Label Content Size 1.6, Show Labels For > Selected clears them, the Shift-drag selects 4 (status bar "Selected: 4") and their four structure labels appear. |
| 04 | Color points by conditions | OK. Conditional gives 4 ranges (18-35.75 ... 71.25-89) and the plot takes them; Bin... 3, Color 2 to grey, patterns <40, 40-65, >65; the grid's AGE cells and the plot's legend follow (26 and 30 blue, 45 to 58 grey). |
| 05 | Compare each group with a control | captions fixed by re-render (0:43, same take). The opening line was a question ("Is each group different from a control group?"), now "Test each group against one control group". "The p-value under the plot" was wrong: p = .0010 is at the top left of the plot; now "The p-value at the top tests all groups at once" (the description too). The maker's question: the two lines pinned at the top cover the Dunnett/Caucasian selectors for ~7 s; left as is, since at the bottom they would cover the Black mark's p-value tooltip, which is what they point to, and the selectors show again at the right-click. Checked: Asian ***, Black ns (adjusted p .4889), Other ***; the table's three rows say the same. |
| 06 | Shade a range of values on a chart | OK. The first region is 36 to 56, which are AGE's first and third quartiles (recomputed: 36, 45, 56; 517 rows); the typed region 65 to 89 shows "65 and over / 85 rows" (recomputed: 85, oldest 89). |
| 07 | Add a note that explains a chart | refilmed (0:44). "Markdown text; ${COLUMN} shows the current row's value" showed over the dialog's default text and was gone a second after the new text arrived; it is now a caption on a hover over the typed text. The maker's question: the note's heading is small (about 15 px at 1080p) in a large empty Markup panel. A "#" heading is no bigger (the default "What's Markdown?" is an H1 of the same size), so it stays; readable at full size. Checked: row 1 (age 26, 74.10 kg) fills the note after OK, the click on row 11's marker gives "age 46, 118.90 kg, 168.974 cm". |

Verdicts: 4 OK (01, 03, 04, 06), 1 captions fixed by re-render (05), 2 refilmed (02, 07).

## Data claims

All true in the frames, and recomputed from data/demo/demog-1000.csv where the clip states a number: 544 of 872 above
BMI 25 (01, the line only claims the BMI); 26 rows in the region (02); 4 selected with labels (03); the three AGE
ranges (04); Asian and Other differ from Caucasian, Black does not (05); quartiles 36 and 56, 85 rows from 65 (06); row
11's values (07).

## Titles and descriptions

Plain, no hype words; each names its menu path and demo dataset. Fixed: 05 (where the p-value is). Checked against
the lead's 23:00 rules: no card, line or description mentions how a video was made (the playlist header of
descriptions.md names the recorder, for maintainers), and the shortest clip is 0:43.

## Systematic

- A line placed before `user enters ... into ...` shows over the field as it was; the typed text arrives at the end of
  the step. A line about typed text belongs on a following hover (07 here, 12-databases/05 the same).
- When a plot is colored by a category, the selection orange can match a category's color (iris versicolor): a clip
  that selects needs another sign of it, such as the status bar's count.
- Lines pinned at the top cover a viewer's own top controls (05); pick the pin by what must stay visible while the line
  shows, which can be a tooltip at the bottom.

## For the lead

- Product/design: the Markup viewer draws Markdown headings (# and ##) at about body size, so a one-line note sits small
  in a large panel (07).
- Product/design: the default selection color is the same orange as the second categorical color, so selected points
  of that category don't change at all (02).
