# Samples changelog

## v.next

* Docker: Raised conda env to Python 3.10 and added security floors (setuptools/wheel/urllib3/pyarrow/jaraco.context/brotli), plus removed stale ensurepip/pkgs-cache copies (VEX)
* Docker: Cleared reported CVEs — added `apt upgrade` for base-image OS packages and upgraded pip/setuptools/wheel in the conda env
* GROK-14287: Dashboards: Fixed errors when opening the "Chemical Space Using tSNE" project
* GROK-14286: Demo Notebooks: Open both the wells table and the notebook view by default

## 1.4.2 (2025-09-22)

* GROK-18923: Northwind: order details by @quantity, @productname, @Country: error
* Samples: Notebooks: Use public python API
* Samples: Deprecated more swager files

## 1.2.0 (2024-10-24)

* Queries auto tests
* Pyodide samples

## 1.1.0 (2024-09-09)

* Dependencies bump

## 1.0.4 (2024-07-18)

* Fixes demo swagger connections, projects, scripts

## 1.0.3 (2024-06-05)

This release focuses on improving feature stability and usability.

* Dependency: datagrok-api >= 1.13.0*
* Fixes demo projects
