# Row Source edge cases — what the features cannot say yet

Written for `row-source-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Questions for owners**

1. A viewer Filter `${AGE} > 44` whose column is then renamed or removed: the formula is treated as empty and the viewer widens to every row without a word (`validateColumnFormulaFields`). Should the formula follow a rename (as formula lines do) and should a removal be reported?
2. Under MouseOverGroup, the viewer whose group is hovered (a pie) closed: keep the group or clear it?

**API or Dart tests to write**

1. Dart regTest: two viewers sharing one Filter formula, one closed while the formula column is computing (`dataFrame.temp` flag).
2. Dart regTest: the Data Manipulation crash test looped over the row sources (it runs on Filtered only).

## Won't do

- Row Source and Filter through a saved layout for the seven viewers: what a layout stores is an API test, and the viewers' own persistence features restore their settings.
