# Markup edge cases — what the features cannot say yet

Written for `markup-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Questions for owners**

1. Html mode content from a shared layout or project with `<img onerror>` / `<script>`: the content is inserted with `NodeTreeSanitizer.trusted`. Is author HTML trusted (d4 rule 8 allows `trusted` for constants only)?

**API or Dart tests to write**

1. Dart regTest: two `{{…}}` expressions on one line (`{{a}} and {{b}}`) are captured as one by the greedy `{{(.*)}}`.

## Won't do

- Long content and numeric bounds: the viewer has no numeric property and long text is plain DOM text.
