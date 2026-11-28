# PMX Studio

PMX Studio helps in vivo pharmacology teams run efficacy studies in mouse tumor models (syngeneic, cell-line
xenograft and PDX) from start to finish: plan a study, follow it while it runs, analyse it when it ends, and compare
it with others. It is made for the scientists who design and read these studies, no statistics background needed: you declare the
model and what you expect each treatment to do, and the app works out the rest.

To open it: **Browse > Apps > Preclinical > PMX Studio**. The Overview shows the three stages (Plan, Monitor,
Analyze), your saved plans and studies, and **Open a study**, which takes a plan, a study or a table of measurements,
picked or dropped anywhere on the page.

## Try it

Each stage's card on the Overview has **Try with demo data** at its foot.

- **Plan**: opens a CT26 study of trametinib, anti-PD-1 and the two combined.
- **Monitor**: the same study at day 9, with animals that need attention planted in it.
- **Analyze**: the finished study, with its mouse PK when the NCA package is installed.
- **Compare**: while you have no studies, **Or compare two demo studies** opens a CT26 and an MC38 study side by side.

The demo files are also in **Browse > Files > Demo > pmx-studio**. All are synthetic except one published
MC38 study (see "Example study" in [DEVELOPER.md](DEVELOPER.md)).

## Plan

- Pick the tumor model, the endpoint volume and the treated groups, each with its kind of therapy and its expected
  effect (TGI or T/C).
- See how many animals each group needs for the power you want, how many to implant, and when the study will end.
- See the power you get with the animals you have, and the smallest effect each group can detect.
- Combinations are recognized from the group names.
- Get the study calendar (to your calendar as `.ics`) and the ethics documents: ARRIVE and
  PREPARE checklists, a study plan and the EU non-technical summary, as Word files.
- A model the atlas has not seen: **+ Add your own model** plans it from your own growth history.
- **File > Save** (Ctrl+S) keeps it in PMX Studio; **File > Export** gives a `.pmxplan` file to send to someone else.

## Monitor

- Follow a running study against its plan: each group's growth and where it is heading.
- See which animals need attention today and why: over or about to pass the humane limit, not measured, losing
  weight, or a value outside its expected range (measure again).
- **Propose groups** randomizes the animals into balanced groups and prints the randomization sheet.
- On a planned review, see what the data say so far.
- Welfare: each group's body weight change and the animals near the limits.

## Analyze

- TGI and T/C per group, with the result in words: slows growth, grows faster, or no clear effect.
- Growth curves, waterfall and spider plots, time to endpoint.
- Each treated animal against its digital twin.
- Combinations: whether the two drugs act more than independently.
- Tolerability: body weight change, weight-loss grades, animals removed early, and each group against the NCI limits.
- With a saved plan: whether the controls grew, the effects held and the study ended as planned.
- Exposure (with the NCA package): each group's effect against its mouse exposure and as a multiple of the patients'.
- Results and protocol text ready to copy, a report, a regulatory (CTD) table and a NONMEM modelling dataset.

## Compare

Tick two studies with measurements or more in **Your studies** on the Overview and press **Compare**, or, in
Analyze, **Study > Compare with other studies…** (the study in Analyze ticked). Back returns to where you came from.

- Studies side by side: growth, dose and effect, exposure and effect.
- Each study's controls against the atlas: did they grow as expected for that model?
- All treated groups in one table.

## Your own data

One row per animal and day, as CSV, TSV or Excel. Column names are recognized in English, with or without a unit
in brackets (`Volume (mm3)`, `Weight [g]`):

| Column | Recognised names (examples) |
|---|---|
| Animal | Animal, Animal ID, Mouse, Subject, ID, Tag, Ear tag |
| Group | Group, Arm, Treatment, Cohort, Dose group |
| Day | Day, Study day, Time, Timepoint, DPI, Week (weeks become days) |
| or Date | Date, Measurement date, Visit date; day 0 is the first measurement, or each animal's Randomization date / Treatment start when given |
| Tumor volume | Volume, Vol, TV, Tumor volume |
| or calipers | Length and Width (and Height); the formula is asked or found from your volumes |
| Body weight (optional) | Weight, BW, Body weight |
| Tumor weight, removal reason (optional, Analyze) | Tumor weight; Removal reason |

Units are converted (cm, cm³, mL, kg, mg); without a unit, mm, mm³ and g are assumed. Wide tables are read too: one
column per day (`Animal, Group, Day 0, Day 3…`), or one column per animal. When a column is ambiguous the app asks;
it never guesses. Before opening, it shows what it read and anything it left out, with the reason.

## Glossary

- **TGI** (tumor growth inhibition): how much less a group's tumors grew than the control's, in percent.
- **T/C**: the treated group's tumor volume over the control's; lower is better.
- **Endpoint**: the tumor volume at which an animal's follow-up ends, set in the plan.
- **Humane limit**: the tumor size (1500 mm³ by default) or weight loss (20%) at which an animal must be removed.
- **Digital twin**: what a treated animal's tumor would have measured untreated, from its starting volume, the
  study's controls and the model's growth.
- **Atlas**: PMX Studio's collection of historical control groups for each tumor model.

## Where the numbers come from

The expected control growth of each model comes from an atlas of historical control groups, mined from the
supplementary data of open-access papers, covering 12 tumor models (7 syngeneic, 3 cell-line xenografts, 2 PDX).
Each treated group is compared with the control, and the results are corrected for comparing several groups at once
(Bonferroni by default; Holm and Benjamini-Hochberg can be chosen). The analysis is the one the plan fixed in advance, so the answer does not depend on
choices made after seeing the data. The source papers and their licences are listed in `files/SOURCES.md`.

## For developers

Functions, file formats, computation engine and how the package is built: [DEVELOPER.md](DEVELOPER.md).
