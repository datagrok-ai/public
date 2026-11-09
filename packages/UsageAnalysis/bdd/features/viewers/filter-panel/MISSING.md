# Filter panel — what the features cannot say yet

## Resolved

- **core#1174, core#1223**: fixed in d4 this round. Removing a column a Combined Boolean flag or an Expression rule uses
  gives the rows back, and an All Columns rule is worked out again over the columns left; a hierarchical card drops a
  removed level, also when a layout names it. Dart tests "d4 | Filters | Removed column".

## Still open

- **public#576**: a click after a custom category order; needs a step that drags an item of the Sort > Custom... list.
- **public#2344**: a hovered category's rows in a scatter plot; needs a `mouse-over rows` count computed from
  `rows.mouseOverRowFunc` when asked.
- **public#605**: the add-filter column picker narrows as you type; needs a step that types into an open picker and
  reads its rows.
- **public#3874**: a multi-value tag set in Column Properties; needs a step for the dialog's tag row and a
  `counts of <caption>` reading.
- **public#2359**: the second card's search icon takes the focus; needs a hand-length press on an icon click
  (`icons.dart` puts the focus back one task later).
- **core#479 (second half)**: the single-value row's close and menu icons; needs a kind for the row (`_singleValueRow`).
- **A range card's min above its max** ("Max should be greater than min"): needs `aria-invalid` and the message on the
  card's fields.
- **core#1224**: Save... under a preset name already taken overwrites it without asking — awaiting a product decision.

## Won't do

- **Edge cases of the 2026-10-26 round dropped in review**: special characters in the card search (Dart "d4 | Filters
  | Categorical search"); is null and is not null (ddt `matcher_test.dart`); a removed column's card going with it
  (`filter-panel-categorical.feature`); a Combined Boolean card on a table without booleans (`combined-boolean-filter.feature`
  claims a boolean column joining later); a closed linked source table (ApiTests "DataFrame: Methods: Link").
- **The free-text input's invalid styling**: the refused rule is claimed by the rules the card lists.
- **A layout naming a removed column on a one-column card**: not reproduced (the other cards come back, no error).
- **Card summary at 19 / 20 categories**: cosmetic.
- **A numeric card of empties, a decimal comma under a comma locale, the text card's `!` and `^`**: not planned.
- **Already covered**: all rows filtered out and a column with one value (`filter-panel-categorical.feature`), Filter out
  missing values (`filter-panel-cloned-view.feature`).
