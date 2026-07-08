# API Tests changelog

## 1.10.3 (WIP)

* Security: rebuilt `apitests-docker-test1` on `python:3.12-alpine` (+ `apk upgrade`, refreshed pip/setuptools/wheel) to clear base-OS CVEs (expat/krb5/openssl/musl) and stale Python tooling.

Add `Dapi: entities.save (polymorphic)` tests covering Project and DataConnection round-trip via `grok.dapi.entities.save`.

Run the Node test runner (`package-test-node`) as ESM via tsx: dynamic `startDatagrok` import, ESM `loadTestFiles`, asset/dayjs shims in `node-test-loader/`, and a shared `test-package` `_package` holder so dapi tests no longer transitively load the browser suite. New `start-node`/`stress-node` scripts.

## 1.10.2 (2026-03-17)

Additonal setOptions tests

## 1.7.17 (2024-08-30)

Add tests for grok.dapi.files.list

## 1.7.9 (WIP)
