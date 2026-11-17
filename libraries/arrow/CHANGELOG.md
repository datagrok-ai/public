# arrow changelog

## v.next

* [core#1407](https://github.com/datagrok-ai/core/issues/1407): Fixed a single-precision float column being built as a double column: the column type follows the array Arrow gives

## 1.1.1 (2026-08-31)

* GROK-20732: Fixed Parquet/Feather nulls reading back as `0` or `1970-01-01`

## 1.1.0 (2026-06-03)

* Fixed `toFeather` column ordering for numeric column names

## 1.0.0 (2026-04-27)

* Extracted from `@datagrok/arrow` with boolean column support
