# DG.Qnum static helpers

**JS API source:** `public/js-api/src/dataframe/qnum.ts:30`

## What we are testing

`DG.Qnum` packs a qualifier (`<`, `=`, `>`) into the two least significant bits of an
IEEE 754 double, so a "qualified number" stays a regular `number` for arithmetic and
sorting. We exercise the helpers that the existing `dataframe.ts` test does not touch:
`Qnum.exact`, `Qnum.less`, `Qnum.greater`, `Qnum.parse`, `Qnum.toString`, and
`Qnum.qualifier`. The first three are convenience wrappers over `Qnum.create` plus a
qualifier constant; the last three round-trip through the Dart side (`grok_Qnum_*`).

## Why it is uncovered

`src/dataframe/dataframe.ts:445` ("qnum" test) only checks
`Qnum.create` / `Qnum.getValue` / `Qnum.getQ`. `src/functions/conversion-functions.ts`
exercises the equivalent **scalar functions** (`ParseQnum`, `Qnum`, `QnumToString`,
`Qualifier`) through the formula engine, but never the static methods on the
`DG.Qnum` class itself. `src/functions/cache.ts` only constructs a value with
`DG.Qnum.less(5)` to feed a cache key.

## Preconditions

- ApiTests package built and published; no DataFrame, server, or filesystem state required.
- The Dart-side bindings `grok_Qnum_Parse`, `grok_Qnum_ToString`, `grok_Qnum_Qualifier`
  are registered (they are — see `core/shared/grok_shared/lib/src/interop/grok_api.dart`).

## Test cases

1. **exact / less / greater set the right qualifier** — for a few sample values, packing
   with each helper preserves the magnitude (to two decimals) and yields the expected
   `Qnum.getQ(...)` constant.
2. **exact / less / greater are equivalent to create+constant** — `Qnum.exact(x)` equals
   `Qnum.create(x, DG.QNUM_EXACT)`, same for less/greater. Confirms no extra logic.
3. **parse round-trip** — `Qnum.parse('10')`, `'<10'`, `'>10'` recover value `10` (to two
   decimals) with qualifiers EXACT / LESS / GREATER. Whitespace tolerance:
   `' < 10'` and `'>  10'` parse the same as their compact forms.
4. **parse rejects garbage** — `Qnum.parse('abc')` and `Qnum.parse('')` are nullish
   (the Dart binding may surface either `null` or `undefined` here, so we check
   `== null` to accept both).
5. **qualifier returns '=', '<', '>'** — `Qnum.qualifier(exact(1.5))` is `'='`,
   `qualifier(less(1.5))` is `'<'`, `qualifier(greater(1.5))` is `'>'`. (The JS-side
   binding maps EXACT to `'='`, not the empty string the raw Dart `QNum.qualifier(int)`
   returns.)
6. **toString shape** — `Qnum.toString(exact(1.5))` has no qualifier prefix; `toString`
   of `less(1.5)` starts with `<`; `toString` of `greater(1.5)` starts with `>`. Each
   strips back to a finite parse via `parseFloat(s.replace(/^[<>=]/, ''))` close to 1.5.

## Out of scope

- The `DG.Column.qnum` constructor and qualified-number aggregation (covered in
  `src/dataframe/dataframe.ts`).
- The `ParseQnum` / `Qnum` / `QnumToString` / `Qualifier` scalar functions exposed via
  the formula engine (covered in `src/functions/conversion-functions.ts` and
  `math-functions.ts`).
- Float-precision behaviour of the two stolen mantissa bits — only validated to two
  decimals, matching the existing `dataframe.ts` style.
