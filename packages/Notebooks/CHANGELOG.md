# Notebooks changelog

## v.next


* GROK-18695: jupyter-notebook container: Moved to python 3.11 and raised security floors (mlflow 3.x, torch 2.13 cpu, tensorflow 2.20+/keras 3, jupyter-server 2.20, nltk 3.9.4, pillow 12.2, urllib3 2.7, and the transitive VEX floors) to clear all fixable CRITICAL/HIGH CVEs
* Moved the Notebooks Playwright E2E suite into the package (playwright/); helpers from @datagrok-libraries/test/src/playwright
* jupyter-notebook container: Migrated to `datagrok-api==0.1.0` (resource-based `api.tables.download`/`upload`); updated "Open as script" parser accordingly

## 1.6.2 (2026-06-16)

* GROK-20225: Fixed HTML preview 404 by exposing only the nginx port (8090) so the proxy reaches `grok_helper`

## 1.6.1 (2026-06-13)

* GROK-19204: Fixed notebook editor losing all cell data when switching tabs
* jupyter-notebook container: Pinned `datagrok-api==0.0.7` to keep notebook boilerplate working

## 1.6.0 (2026-06-12)

* jupyter-notebook container: Merged the datagrok/jupyter_notebook base image build into the package Dockerfile, making it self-contained (no dependency on the separately built base image)

## 1.4.0 (2025-07-28)

Bug fixes, datagarok-api >= 1.26.0*

## 1.3.0 (2025-04-08)

Bug fixes, Docker container with jupyter_notebook server

## 1.2.0 (2024-10-28)

Fixed corrupted styles for the platform that used to appear after opening a notebook.

## 1.0.1 (2023-07-24)

Initial release of notebooks package

*Dependency: datagarok-api >= 1.0.0*

### Features

Package implements the [Jupyter Notebooks](https://jupyter.org/) integration.
