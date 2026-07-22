# DB Tests changelog

## v.next

* Tests: annotated API-only categories (queries, TableQueryBuilder, server cache, benchmarks, DB annotations, provider connectivity) with `node: true` — `grok test` now runs them headless in Node; the data-sync view test and the IndexedDB-backed client cache stay in the browser
* Tests: the query `-- test:` auto tests (~200) also run headless now that the Node runtime provides the shared `OpenFile` fallback — the whole package needs the browser for just 5 tests

* Tests: grouped connectivity tests by provider (`Providers: <dataSource>`) via `initPackageTests`; removed the dead `tests/categories.ts` whose provider categories held a `before` but no tests.
* Tests: skip all ClickHouse tests — the connectivity check (`initPackageTests` skipReason) and the 31 query round-trip tests (`skip:` on the `-- test:` annotations in `clickhouse-*.sql`). The ClickHouse demo DB is down (failing on dev + CI since ~2026-07-12); remove the `skip:` tokens to re-enable once it's restored.
* Security: pinned the DB test image to `postgres:17-bookworm` and added `apt-get upgrade` to clear stale Debian base CVEs (gnutls28/perl/glibc).

## 1.3.0 (2025-07-28)

* Datagrok api dependency >= 1.26.0

## 1.1.3 (2025-02-18)

* Test fixes

## 1.1.1 (2024-10-21)

* Test fixes
* Eslint version update
* Visible only for developers

## 1.1.0 (2024-09-09)

* Datagrok api dependency > 1.21.1

## 1.0.14 (2024-08-01)

Fixes flapping tests
Implicitly marks benchmarks

## 1.0.13 (2024-05-31)

Fixes connections tests initialization
Fixes data sync test

## 1.0.12 (2024-05-22)

Fixes in server cache tests

## 1.0.8

Initial version of DbTests package

*Dependency: datagarok-api >= 1.16.0*

### Features

DB Connectivity test suite
