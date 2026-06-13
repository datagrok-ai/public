# Notebooks changelog

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
