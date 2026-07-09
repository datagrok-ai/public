# arrow changelog

## v.next

* GROK-20379: Added `fromFeather` decoding for 128/256-bit decimals
* GROK-20379: Added nested types (list, struct, map) to `fromFeather` as Java-formatted strings
* GROK-20379: Added the `narrowFloatsToFloat32` option to `fromFeather`
* GROK-20379: Fixed 64-bit ints mapping to `int` when values happened to fit

## 1.1.0 (2026-06-03)

* Fixed `toFeather` column ordering for numeric column names

## 1.0.0 (2026-04-27)

* Extracted from `@datagrok/arrow` with boolean column support
