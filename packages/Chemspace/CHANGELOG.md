# Chemspace changelog

## v.next

* App: Surfaced search errors instead of silently showing an empty table, and guarded against stale responses from rapid input changes
* Prices panel: Fixed caching — country switches now reuse rendered results instead of refetching
* Get Chemspace Prices: Offers with missing USD pricing now fall through to the next offer instead of being stored as undefined
* queryMultipart: URL-encoded query params (raw joining broke on spaces and commas)
* queryMultipart: Rewrote retry loop — a successful 401 refresh no longer throws, and refreshed tokens are actually used on retry

## 1.2.3 (2025-05-06)

Fixed token request

## 1.2.2 (2025-04-22)

Removed swagger file

## 1.2.1 (2025-03-29)

Moved application under Chem section in browse tree

## 1.2.0 (2025-03-29)

Datagrok api >= 1.25.0

## 1.1.6 (2025-02-25)

Updated README

## 1.1.5 (2025-02-17)

### Features

* Updated tests

## 1.1.4 (2024-11-18)

### Features

* Added function to add Chemspace ids via 'Add new column' dialog
* Highligh chemspace_id semantic type

### Bug Fixes

* Improved errors reporting

## 1.1.3 (2024-11-18)

### Bug Fixes

* GROK-17022: Chemspace: opening Similar in context panel causes errors

## 1.1.2 (2024-10-16)

* Fixed cache invalidating

## 1.1.1 (2024-08-28)

* Caching results

## 1.1.0 (2024-07-31)

* Updated to Chemspace API version 4.0

## 1.0.5 (2024-07-24)

* Dependency: datgarok-api >= 1.20.0

## 1.0.4 (2024-05-23)

### Bug fixes

* package complience to contemporary platform version

## 1.0.3 (2023-12-19)

### Bug fixes

* Fixed searching for molblocks

## 1.0.2 (2023-07-24)

This release focuses on improving package usability.

*Dependency: datgarok-api >= 1.15.2*

### Features

* Added similarity score for compounds in the similarity panel.
