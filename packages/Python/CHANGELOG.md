# Python changelog

## v.next

* Docker (datagrok/python base): Cleared reported CVEs — added `apt upgrade` for OS packages and upgraded pip/setuptools/wheel; also benefits images built FROM this base (e.g. Bio, Notebooks)
* Docker (datagrok/python base): Synced the stale ensurepip bundled pip wheel with the installed pip (VEX: clears pip 23.0.1 CVE rows while keeping `python -m venv` functional)

## 0.0.1 (2024-11-15)