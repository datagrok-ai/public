# Example study: MC38, LD4172 and anti-PD-1

`mc38-ld4172-anti-pd1.csv` is a published mouse study, offered in the empty Monitor and Analyze as a table to try.

## Source

Yu X. et al. (13 authors; last author Wang J.), "Development of a RIPK1 degrader to enhance antitumor immunity",
Nature Communications 15 (2024), https://doi.org/10.1038/s41467-024-55006-2, Supplementary Fig. 4A–B, Source Data
file (`41467_2024_55006_MOESM4_ESM.xlsx`, sheet "Supplementary Figure4"). Source data are also available in
Figshare, https://doi.org/10.6084/m9.figshare.27629700.

Licence: Creative Commons Attribution 4.0 International (CC BY 4.0),
https://creativecommons.org/licenses/by/4.0/.

The study: MC38 tumors in mice, four groups of 8, randomized and treated from day 7 after inoculation: vehicle;
LD4172 20 mg/kg i.p. daily; anti-PD-1 100 µg i.p. every 3 days; and the combination of both. The Methods state tumors
of ~100 mm³ at that point; the file's volumes on that day: median 64, range 38–103 mm³. Tumor volume (mm³) and body
weight (g) on days 7, 10, 13, 16 and 19 after inoculation: days 0, 3, 6, 9 and 12 in the file.

This is the one study the package distributes animal by animal, under its licence. Its control group is also among
the MC38 studies the package's growth model was built from (`files/SOURCES.md`): for the control band it is not a new
study.

## Changes made

- Reshaped to long format: one row per animal and day (`animal`, `group`, `day`, `volume`, `body_weight`).
- Day re-expressed as days from randomization: article day 7 = day 0.
- Animal labels are ours (V1–V8, P1–P8, L1–L8, C1–C8), assigned in the order of the source columns.
- Group names expanded with the dose, route and schedule given in the article's Methods (µg written "ug").
- Volume and body weight joined by column position: the article gives no animal identifiers; each weight is paired
  with the tumour volume in the same column position of its group.
- Volumes and weights unchanged.
