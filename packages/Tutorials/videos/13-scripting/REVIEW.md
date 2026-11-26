# Review: 13-scripting (Scripting and automation)

Fresh-eyes review (2026-11-22, ~23:15-23:20) of the 6 delivered clips: feature, descriptions.md block, and contact
sheets at the key lines and results of each delivered take (each matched its latest take's length; 04 after its
maker's rework under the lead's newest rules). Judged on: is the one idea
clear to someone new to Datagrok, is each line true when it shows, plain wording, jumps, nothing about how the video
was made, and no clip under ~35 s.

| NN | Clip | Verdict |
|---|---|---|
| 01 | The console | OK. The Add New Column dialog's BMI column is logged as `AddNewColumn("demog", "${WEIGHT} / Pow(${HEIGHT} / 100, 2)", "BMI")`; the pasted command with `${WEIGHT} * 2.2` adds "Weight in lb" (row 1: 161.04 = 73.20 × 2.2); "SelectRandomRows" alone opens its form, Count 100, and the console shows `selectRandomRows("demog", 100, 42, true)` with "Selected: 100" in the status bar. |
| 02 | Write a JavaScript script | OK. The header lines are typed key by key, so each line shows while its text appears; Run's form has the two captioned inputs, Height 180 gives bmi 21.60 (70 / 1.8²) under the code; SAVE, and the Scripts search finds the new card. The Scripts gallery shows test-looking scripts for ~4 s (stand clutter, the maker's note). |
| 03 | Run a script on your table | OK. The header (a dataframe input, two numerical column inputs) is pointed at with its line; the run form fills Height and Weight by name once demog is picked, and the clip picks them again; demog gets a BMI column (row 1: 28.42 = 73.2 / 1.60484²). Same gallery clutter for a few seconds. |
| 04 | Run a saved script without opening it | OK (the maker's rework, 0:41). The card's context panel: Details (BMI from height and weight, inputs height, weight, output bmi), the Script pane with the code, the Run pane's form (175/70 by default) set to 160 and 55, RUN, and the Variables panel reads bmi 21.48 (55 / 1.6²). The feature file is still named 04-run-a-saved-script-again.feature while the title says "without opening it" (maintainer-facing only; renaming needs a file removal, left for the morning). |
| 05 | Learn from the sample scripts | OK. JS API Samples > grid > color-coding opens as a script; the height line is changed to '20-160' green and '160-190' dark; F5 opens "demog 10000" where heights under 160 are green, 160 to 190 dark and the few over 190 uncolored, as the changed line says. |
| 06 | Give a script a friendly form | OK. Each header option has its line while pointed at; the run form shows the column picker, the 1-to-4 slider set to 2.5 and the Then list; 114 rows selected, recomputed from data/demo/demog.csv (WEIGHT mean 78.74, sample SD 19.92, 114 values more than 2.5 SD away). |

Verdicts: 6 OK, no changes by the reviewer.

## Data claims

All true in the frames: 161.04 (01), 21.60 (02), 28.42 (03), 21.48 (04), the green/dark/uncolored heights (05), 114
rows (06).

## The lead's newest rules

No card sentence, line or description mentions how the video was made; the shortest clip is 04 at 0:41.

## Systematic

- The Scripts gallery (Browse > Platform > Functions > Scripts) lists the stand's test scripts (StressMutateStep,
  MockPipeline3, HiddenOutputsTest, ViewersHookModel ...) for a few seconds before each clip's search (02, 03, 04, 06).
- A saved button's tooltip ("Save script (Ctrl+S)") stays up after the pointer leaves for the next press (02).

## For the lead

- Stand clutter: the test scripts in the Scripts gallery (above), the same cleanup question as the dashboards gallery.
