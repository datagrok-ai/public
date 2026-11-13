# Pie chart — what the features cannot say yet

## Won't do

- Max Radius and Start Angle typed past their bounds: core#1200's clamp tests (`../MISSING.md`).
- A Filter formula no row passes, and a click on the empty pie: no rows is the Dart Data Manipulation loop; a click that
  changes nothing is no claim.
- A renamed or removed category column, no rows, deleted rows, a NaN or ±Infinity angle column, a one-row table,
  200000 categories: the rules every viewer shares (`../MISSING.md`) and xamgle `d4 | Viewers | Non-finite values | Pie
  chart`.
