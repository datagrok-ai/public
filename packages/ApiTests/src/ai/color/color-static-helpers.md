# DG.Color static helpers

**JS API source:** `public/js-api/src/color.ts:30`

## What we are testing
Pure JS static helpers on `DG.Color` that pack and unpack ARGB-formatted
integers without crossing the Dart interop boundary: channel extractors
(`a`, `r`, `g`, `b`), the inverse packer `argb`, `setAlpha`, the HTML
round-trip pair `fromHtml` / `toHtml`, the `rgb(...)` formatter `toRgb`, the
linear `scale` mapper, the `hexToPercentRgb` decomposer used by the RDKit
substruct highlight, and the categorical-palette wrap-around accessor
`getCategoricalColor`.

## Why it is uncovered
`src/utils/color.ts` exercises a sibling utility module
(`@datagrok-libraries/utils/src/color`) and only touches one method on
`DG.Color` — `DG.Color.toHtml`. None of the channel extractors,
`argb` / `setAlpha`, `fromHtml`, `toRgb`, `scale`, `hexToPercentRgb`, or
`getCategoricalColor` are referenced anywhere under
`public/packages/ApiTests/src/` outside `src/ai/`.

## Preconditions
- None. All helpers under test are pure static functions on `DG.Color`
  except `getCategoricalColor`, which reads the categorical palette via the
  Dart interop layer — the test asserts the wrap-around modulo behaviour
  against the palette returned by `DG.Color.categoricalPalette`, so any
  palette length is tolerated.

## Test cases
1. **channel extractors** — for the ARGB integer `0x12345678`, expect
   `Color.a` == `0x12`, `Color.r` == `0x34`, `Color.g` == `0x56`,
   `Color.b` == `0x78`.
2. **argb round-trip** — `Color.argb(0x12, 0x34, 0x56, 0x78)` reproduces
   `0x12345678`, and feeding the result back through the channel extractors
   returns the original components.
3. **setAlpha preserves rgb** — `Color.setAlpha(0xFF112233, 0x80)` equals
   `0x80112233`; the original color is not mutated.
4. **fromHtml / toHtml round-trip** — `Color.fromHtml('#00bfff')` produces
   the integer whose lower 24 bits are `0x00bfff`, and feeding that integer
   into `Color.toHtml` yields `'#00bfff'`.
5. **toRgb formatting** — `Color.toRgb(0xFF00BFFF)` equals
   `'rgb(0,191,255)'`.
6. **scale linear mapping** — `Color.scale(5, 0, 10)` is `0.5`,
   `Color.scale(0, 0, 10)` is `0`, `Color.scale(10, 0, 10)` is `1`, and a
   degenerate range (`min == max`) returns `min` instead of dividing by
   zero.
7. **hexToPercentRgb happy path** — `'#00bfff'` decomposes to
   `[0/256, 0xbf/256, 0xff/256, 0.3]` (the default alpha used by the RDKit
   highlight).
8. **hexToPercentRgb with explicit alpha** — `'#00bfff80'` decomposes to
   `[0/256, 0xbf/256, 0xff/256, 0x80/256]`.
9. **hexToPercentRgb invalid input** — `'not-a-color'` returns `null`.
10. **getCategoricalColor wrap-around** — for any non-empty palette of
    length `N`, `getCategoricalColor(i)` equals
    `categoricalPalette[i % N]` for `i = 0`, `i = N - 1`, `i = N`, and
    `i = 2*N + 3`.

Negative cases beyond `hexToPercentRgb('not-a-color')` are skipped: the
remaining helpers are pure arithmetic on `number` inputs with no defined
failure mode (passing `NaN` or out-of-range integers is observable but not
documented and would lock in undefined behaviour).

## Out of scope
- Dart-backed helpers (`getCellColor`, `getCategoryColor`,
  `getContrastColor`, `getRowColor`, `scaleColor`, `highlight`, `darken`)
  — these need real cells/columns and belong in a follow-up scenario under
  the same `src/ai/color/` folder.
- The named-color constants (`Color.gray`, `Color.blue`, etc.) — trivial,
  not worth a test each; would belong in a snapshot-style scenario if
  anything.
- Color coding on grid columns — already covered by
  `src/grid/color-coding.ts`.
