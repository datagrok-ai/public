# Docking changelog

## v.next

* Tests: The AutoDock tests are skipped, with the reason, when the server reports Grok Spawner switched off or unreachable instead of each waiting 5 minutes for a container that never starts
* Docker: autodock — raised build-only base env floors (urllib3>=2.8, msgpack>=1.2.3, py-rattler>=0.28.1 for the zip crate CVE) and removed the base env's pip with its vendored urllib3/msgpack/setuptools (VEX)

## 1.3.4 (2026-11-03)

* Docker: rebuilt the autodock image without the install-only pip/setuptools/wheel/certifi in the Python 2.7 env (VEX)

## 1.3.3 (2026-07-12)

* Docker: autodock — raised build-only base env floors (urllib3>=2.7, msgpack>=1.2.1) and removed stale py2.7 ensurepip bootstrap wheels (VEX)

## 1.3.2 (2025-04-24)

* Docking: Demo improvements

## 1.1.2 (2024-10-21)

### Bug Fixes

* Docking: Increased timeout for clinfo test

## 1.1.1 (2024-10-18)

### Bug Fixes

* Docking: Fixed dockLigandColumn test

## 1.0.4 (2024-08-23)

### Features

* Docking: Add demo
* Docking: Ensure consistent naming and ordering of parameters

## 1.0.3 (2024-04-15)

Add description to the Docking demo card.

## 1.0.2 (2024-04-12)

### Bug Fixes

* Upgraded BiostructureViewer package to a newer version

## 1.0.1 (2024-04-08)

### Bug Fixes

* [#2741](https://github.com/datagrok-ai/public/issues/2741): Docking: Show placeholder for failed docking run

## 1.0.0 (2024-03-19)

### Features

* Split from BiostructureViewer package
* [#2491](https://github.com/datagrok-ai/public/issues/2491): Docking: Autodock-GPU integration
* [#2568](https://github.com/datagrok-ai/public/issues/2568): Docking: Input UI
* [#2732](https://github.com/datagrok-ai/public/issues/2732): Docking: Handle docking error for a particular ligand

### Bug Fixes

* Fix AutoDock for pose_count param
* Fix AutoDock convert, add tests
* Fix AutoDock add clinfo test
