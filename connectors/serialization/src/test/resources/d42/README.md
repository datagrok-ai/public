# d42 cross-language fixtures (GROK-20345 / WO-2, GROK-20346 / WO-3)

`*.d42` are **Dart-written** d42 blobs; each `*.expected.json` sidecar records, per
column, the wire type, the on-wire encoder id, and the expected values. They prove that
the Java d42 reader (`serialization.DataFrame.fromByteArray`) decodes Dart output
value-for-value, including best-encoder selections the Java writer never emits.

`D42DartFixtureTest` decodes every `.d42` here and asserts each cell against its sidecar
(float64/qnum bit-exact, float32 at float precision), the observed encoder id against the
sidecar (so a ddt cost-model drift fails loudly), and — in `testIdInventoryCovered` — that
the **full** ARCHITECTURE §1.5 encoder-id matrix is exercised (`int 1/2/3/4`,
`string 0/1/2/3`, `bool 1`, `datetime 1/2/3`, `double 1/2/3/4/5`, `bigint 1/2/3`, `qnum 1`,
`byte_array 1`, `dataframe 1`). `testUnknownEncoderIdThrows` pins the loud-error contract
against the hand-corrupted `corrupt_unknown_encoder.d42` (no sidecar).

## Fixture families

- **forced_*** — each column's encoder id is pinned directly (the generator writes the id
  and calls that encoder), so every id has a deterministic fixture regardless of the cost
  model. Covers the encoders the natural cost model would never pick (fcp, legacy
  rle/raw floats, string prefixes/squash/zlib, datetime raw/int component forms, bigint
  capped, qnum, byte_array, recursive dataframe).
- **natural_*** — real best-encoder selection (`col.meta.encode`), engineered so the cost
  model picks the target encoder (arithmetic → `int:pattern`; low-cardinality →
  `int:bitIntList`; 3-category → `string:categories` with `bitIntList` indices; distinct
  doubles → `float:raw64`; a plain bigint column → `bigInt:list` (id 2); a production-shaped
  mixed frame; a 0-row frame).
- **corrupt_*** — hand-corrupted (a valid int frame with the encoder id overwritten by an
  unregistered value); excluded from the parameterized decode + inventory, asserted only by
  `testUnknownEncoderIdThrows`.

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
