# Markup edge cases — what the features cannot say yet

Written for `markup-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Questions for owners**

None. Html mode content from a shared layout or project runs its scripts (`NodeTreeSanitizer.trusted`): core#1239, no scenario (a script run is not observable from a step).

**API or Dart tests to write**

1. Dart regTest: two `{{…}}` expressions on one line (`{{a}} and {{b}}`) are captured as one by the greedy `{{(.*)}}`.

## Won't do

- Long content and numeric bounds: the viewer has no numeric property and long text is plain DOM text.
