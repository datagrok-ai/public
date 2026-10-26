# Pie chart edge cases — what the features cannot say yet

Written for `pie-chart-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Questions for owners**

1. The category column removed from the table while the pie shows it (and renamed): fall back to another categorical column, or an empty state with a message?
2. A layout whose category or angle column the table lacks: the pivot reports the missing column; what should the pie do?
3. A Segment Angle Column holding NaN or ±Infinity: refuse with a message like the negative and all-zero cases?
4. Max Radius typed below 50 or above 500, or cleared; Start Angle −30 or 400: clamp, refuse, or auto?

**API or Dart tests to write**

1. Dart regTest: more than 200000 categories shows "too many categories" (`pie_chart_core.dart`, `validateState`).

## Won't do

- The other core#982 fixes (On Click switched to Select and back, the category changed while the pie filters, an unsupported Category Map from a layout, the table switch keeping the columns): each has a Dart regTest in `viewers_test.dart` ("Fixed issues | Pie chart …"); the feature keeps the one whose gesture is a wedge click under a Filter formula.
