# CVM Tests changelog

## v.next

* Security: rebuilt `cvmtests-docker-test1` (`python:3.12-alpine`) and `cvmtests-docker-test2` (`python:3.11-alpine`) on current bases (+ `apk upgrade`, refreshed pip/setuptools/wheel) to clear base-OS CVEs (expat/krb5/openssl/musl) and stale Python tooling.
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
