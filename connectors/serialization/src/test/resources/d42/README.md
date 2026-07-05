# d42 cross-language fixtures (GROK-20345 / WO-2)

`*.d42` are **Dart-written** d42 blobs; each `*.expected.json` sidecar records, per
column, the wire type, the on-wire encoder id, and the expected values. They prove that
the Java d42 reader (`serialization.DataFrame.fromByteArray`, WO-1) decodes Dart output
value-for-value, including best-encoder selections the Java writer never emits.

`D42DartFixtureTest` decodes every `.d42` here and asserts each cell against its sidecar
(float64 bit-exact, float32 at float precision), the observed encoder id against the
sidecar (so a ddt cost-model drift fails loudly), and — in `testIdInventoryCovered` — that
the WO-1 decoder set (`int 1/2/3/4`, `string 0`, `bool 1`, `datetime 3`, `double 1/5`,
`bigint 1`) is exercised.

## Fixture families

- **forced_*** — each column's encoder id is pinned directly (the generator writes the id
  and calls that encoder), so every WO-1 id has a deterministic fixture regardless of the
  cost model. `bigint` (`bigInt:raw`, id 1) and `float:raw`/`float:raw64` (ids 1/5) live
  only here — the natural cost model would pick ids WO-1 does not decode yet (WO-3).
- **natural_*** — real best-encoder selection (`col.meta.encode`), engineered so the cost
  model picks the target encoder (arithmetic → `int:pattern`; low-cardinality →
  `int:bitIntList`; 3-category → `string:categories` with `bitIntList` indices; distinct
  doubles → `float:raw64`; a production-shaped mixed frame; a 0-row frame). Deliberately
  no `bigint` (its natural encoder is `bigInt:list`, id 2 — WO-3).

## Regeneration (byte-identical-copy rule)

These files are **byte-identical copies of the generator output** — do not hand-edit. To
regenerate after a legitimate format/fixture change:

```
cd core/shared/ddt
pub run test test/serialization/d42_fixture_generator_test.dart
```

The generator (`core/shared/ddt/test/serialization/d42_fixture_generator_test.dart`) writes
the `.d42` + `.expected.json` pairs straight into this directory, round-trips each fixture
through the Dart reader (guarding against generator bugs), and asserts the natural family
picks the intended encoder. Serialization is deterministic (metadata ids are pinned and the
blob tail carries a minimal fixed JSON), so an unchanged format leaves the tree clean. A
diff here after regeneration is the loud signal that the format or the ddt cost model moved
— review it, then re-run `D42DartFixtureTest`.

Fixtures are kept under 100 KB each (they live in git).
