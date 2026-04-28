# DG.StringUtils static helpers

**JS API source:** `public/js-api/src/helpers.ts:185`

## What we are testing

`DG.StringUtils` exposes two static helpers used across the platform: `hashCode(s)` — a
pure-JS 32-bit string hash with the classic `(h<<5) - h + chr` recurrence — and
`camelCaseToSentence(s, opts?)` — a Dart-backed converter that splits a camelCase
identifier into a space-separated phrase, with toggles for capitalising the first word,
each subsequent word, and conjunctions (`and`, `or`, `than`, `if`, `but`, `so`, `as`,
`that`). Both are deterministic and have no server dependency beyond the in-page Dart
runtime.

## Why it is uncovered

Grepped `public/packages/ApiTests/src/` (excluding `src/ai/`) for `StringUtils`,
`hashCode`, and `camelCaseToSentence` — zero hits. No ApiTests file exercises this class
today.

## Preconditions

- ApiTests package loaded; `DG.StringUtils` resolves (i.e. js-api bundle is the running one).
- No DataFrame or server state needed.

## Test cases

1. **hashCode empty string** — `DG.StringUtils.hashCode('')` returns `0`.
2. **hashCode single char** — `DG.StringUtils.hashCode('a')` returns `97`
   (`(0<<5) - 0 + 97`).
3. **hashCode is deterministic and discriminates inputs** — same input yields the same
   number across calls; a small fixed set of distinct strings yields distinct hashes.
4. **hashCode returns a 32-bit signed integer** — for a long input, the result is an
   integer in `[-2^31, 2^31)` (verifies the `|0` truncation).
5. **camelCaseToSentence default** — `'helloWorld'` becomes `'Hello World'` (first
   word capitalised; subsequent words preserved as `splitCamelCase` yields them).
6. **camelCaseToSentence capitalizeFirst=false** — `'helloWorld'` with
   `{capitalizeFirst: false}` becomes `'hello World'`.
7. **camelCaseToSentence capitalizeConjunctions** — `'fooAndBar'` becomes `'Foo and Bar'`
   by default (the conjunction `and` is lower-cased) and `'Foo And Bar'` with
   `{capitalizeConjunctions: true}`.
8. **camelCaseToSentence pass-through** — inputs that are fully upper-case (`'ID'`) or
   already contain a space (`'hello world'`) are returned unchanged.

Negative cases: the API has no defined failure mode for these helpers (no thrown
exceptions, no nullable returns documented), so coverage is limited to observable
transformations.

## Out of scope

- `FormulaLinesHelper` and `AnnotationRegionsHelper` (separate scenarios; they touch
  viewer/dataframe metadata and deserve their own folder).
- The `capitalizeNext` toggle — for typical camelCase input the post-split second word
  is already capitalised, so the flag is observationally a no-op; exercising it
  meaningfully needs non-camelCase input and a separate scenario.
- Round-trip with non-ASCII identifiers — the underlying `splitCamelCase` uses ASCII
  A-Z bounds, so behaviour is undefined for unicode and we do not assert on it.
