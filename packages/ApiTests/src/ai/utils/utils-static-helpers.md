# DG.Utils static helpers

**JS API source:** `public/js-api/src/utils.ts:171`

## What we are testing

Pure-JS static helpers on `DG.Utils`: `firstOrNull`, `identity`, `replaceAll`,
`isEmpty`, `nullIfEmpty`, `getJsonValueType`, `jsonToColumns`, `randomString`.
These are deterministic helpers that do not round-trip through Dart, so they
can be exercised without a server. `random<T>(items)` is intentionally skipped
because it returns a non-deterministic element.

## Why it is uncovered

`src/utils/progressIndicator.ts` covers `ProgressIndicator`, not the `Utils`
class. `DG.Utils.randomString(...)` is *used* as a fixture in `src/dapi/*.ts`
but is never asserted. A grep of `src/` (excluding `src/ai/`) for
`DG.Utils.firstOrNull|identity|replaceAll|isEmpty|nullIfEmpty|getJsonValueType|jsonToColumns`
returns no hits. The sister area `src/ai/utils/string-utils-helpers.{md,ts}`
(round 3) covers `DG.StringUtils` only, so the `Utils` class itself is
genuinely uncovered.

## Preconditions

- ApiTests package built and published to `dev`.
- No data, server state, or external packages required.

## Test cases

1. **firstOrNull-on-empty-iterable** — pass `[]`, expect `null`.
2. **firstOrNull-on-array** — pass `[10, 20, 30]`, expect `10`.
3. **firstOrNull-on-set** — pass `new Set(['a', 'b'])`, expect `'a'` (insertion
   order is well-defined for `Set`).
4. **identity-zero-and-positive** — `identity(0)` is an empty `Uint32Array`;
   `identity(5)` equals `[0, 1, 2, 3, 4]` and is a `Uint32Array`.
5. **replaceAll-replaces-every-occurrence** — `'a-b-a-c-a'` with `'a' -> 'X'`
   gives `'X-b-X-c-X'`. Empty `replace` strips the search token. Search token
   absent from input → string returned unchanged.
6. **isEmpty-handles-null-undefined-empty** — `isEmpty(undefined)`,
   `isEmpty(null as any)`, `isEmpty('')` all `true`; `isEmpty(' ')` and
   `isEmpty('x')` are `false`.
7. **nullIfEmpty-returns-null-or-original** — `nullIfEmpty('')`, `nullIfEmpty(undefined)`
   are `null`; `nullIfEmpty(' ')` returns `' '` (whitespace is not empty).
8. **getJsonValueType-known-and-falsy-and-throw** — `'hi' -> 'string'`,
   `42 -> 'double'`, `true -> 'bool'`. Falsy inputs (`null`, `undefined`,
   `''`, `0`, `false`) all return `null` because of the `!x` short-circuit.
   An object input throws.
9. **jsonToColumns-expands-string-column** — DataFrame with one string column
   `j` of `{"a":1,"b":"x"}` and `{"a":2,"b":"y"}` rows; after
   `Utils.jsonToColumns(col)`, the frame has columns `a` (double) and `b`
   (string) with the expanded values.
10. **randomString-length-and-alphabet** — `randomString(0)` is `''`;
    `randomString(20)` has length 20 and contains only `[A-Za-z0-9]`; two
    calls with the same length are not guaranteed equal but must respect the
    alphabet and length contract.

## Out of scope

- `Utils.random(items)` — non-deterministic, only worth a smoke check; skip.
- `Utils.openFile`, `openFileBytes`, `download` — DOM/browser dialogs,
  belong in widget/file-viewer integration tests.
- `Utils.loadJsCss`, `Utils.streamToObservable`, `Utils.detectColumnHierarchy`,
  `Utils.executeTests`, `Utils.createResultsCsv` — Dart/server/test-harness
  integration; out of scope for a static-helper sweep.
