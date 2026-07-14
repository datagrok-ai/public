# CVM Tests changelog

## v.next

* Tests: fixed Docker tests — the container-name filter never matched. Platform registers package containers as `kebab(package.name)-<dockerfileFolder>` (`cvm-tests-cvmtests-docker-test1/2`); the test queried `cvmtests-Cvmtests-...`, so `before()` got `undefined` and all 5 Docker tests failed. Also added a clear not-found error instead of a cryptic undefined cascade.
* Tests: fixed `Column list` script tests — JS `column_list` input is a name array (`string[]`), so use `cols[0]` not `cols.toList()[0]`; Grok script used an undefined `columns` var and numeric indexing (Grok lists expose `first`/`last`/`length`, not `[i]`), now `cols.first`.
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
