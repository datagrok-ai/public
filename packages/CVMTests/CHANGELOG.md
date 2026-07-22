# CVM Tests changelog

## v.next

* Tests: Added `Celery: node worker` suite — JS queue functions (`meta.queue`) executed server-side in the celery Node worker: scalar/string/dataframe round-trips, error propagation, progress, cancellation, dapi token pass-through, and a custom `dockerfiles/node-worker` container

* Tests: env scripts now declare `#meta.timeout: 600` (the 300s server exec default fired at ~318s under merged-stage load) and env-test client budgets sit above it at 660s — real hangs still surface, a busy stand doesn't flake
* GROK-20445: Env tests: dropped the obsolete `numpy < 2` assertion (base env now resolves numpy-2-native pyarrow); made the two DG-idiom nodejs scripts (Column list, Lines count) runtime-agnostic until the Node js-api runtime lands on master
* Scripts: Migrated nodejs test scripts to the js-api (`grok`/`DG` globals); dataframe-js/axios are no longer injected into nodejs scripts (breaking — legacy scripts must `require()` them explicitly or move to `DG.DataFrame`)

* Tests: fixed Docker tests — the container-name filter never matched. Platform registers package containers as `kebab(package.name)-<dockerfileFolder>` (`cvm-tests-cvmtests-docker-test1/2`); the test queried `cvmtests-Cvmtests-...`, so `before()` got `undefined` and all 5 Docker tests failed. Also added a clear not-found error instead of a cryptic undefined cascade.
* Tests: fixed `Column list` script tests — JS `column_list` input is a name array (`string[]`), so use `cols[0]` not `cols.toList()[0]`; Grok script now `DeleteColumns(df, Named([cols.first]))` — a runtime list value isn't coerced to a column filter (only parse-time list literals are), so the prior `DeleteColumns(df, [cols.first])` threw `Class 'List' has no instance method 'makePredicate'`. Wrapping in `Named(...)` builds the predicate explicitly.
* Security: rebuilt `cvmtests-docker-test1` (`python:3.12-alpine`) and `cvmtests-docker-test2` (`python:3.11-alpine`) on current bases (+ `apk upgrade`, refreshed pip/setuptools/wheel) to clear base-OS CVEs (expat/krb5/openssl/musl) and stale Python tooling.
* Docker: cvmtests-docker-test2 — raised Quart (>=0.20) and Werkzeug (>=3.1.6) floors to clear their CVEs (VEX)
* Added datagrok-celery-task integration tests via the python/ celery worker

# 1.4.0 (28-07-2025)

* Dependency: datagarok-api >= 1.26.0*

# 1.2.0

* Test fixes
* Dependency: datagarok-api >= 1.24.0*

# 1.1.1

* Test fixes

# 1.1.0

* Version bump

# 1.0.10 (WIP)

* Skip some tests, dependencies bump

# 1.0.9

* Fixed scripting tests. Added test for int columns

## 1.0.6 (WIP)
