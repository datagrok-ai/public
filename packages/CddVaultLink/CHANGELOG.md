# CDD Vault Link changelog

## v.next

* Tab loading: Unified all tabs on a preview-then-"Load all" pattern. Initial open always shows the first 100 rows via the sync endpoint; if more exist, a ribbon appears with a "Load all" button that runs the async endpoint in the background and swaps the dataframe on completion. Removed the old `createCDDTableViewWithPreview` helper.
* API: Added sync `queryCollections` and `queryBatches` wrappers to pair with their async counterparts.
* Types: Added `molecule` (parent `Molecule`) to `Batch`; added `registration_type`, `registration_form`, `num_aromatic_rings`, `p_k_a_basic`, `bbb2_score` to `Molecule`; made `Molecule.batches` optional (absent when a Molecule is nested inside a Batch).
* Renamed search functions: `cDDVaultSearch2` → `cDDVaultSearch` (now the canonical sync counterpart to `cDDVaultSearchAsync`, with a `page_size` parameter); the old 25-parameter `cDDVaultSearch` moved to `cDDVaultSearch2`.

## 1.1.0 (2026-03-16)

### Bug fixes

* GROK-19523: AppTreeBrowser decorator: Fixed behavior, so it finds the app properly, fixed usages

## 1.0.5 (2025-12-24)

### Bug fixes

* GROK-18964 CDD Vault: Cannot get statistics for vault DataGrok test vault: Failed to fetch

## 1.0.4 (2025-05-28)

* Some UI fixes
* Having only one view opened in a preview at once

## 1.0.3 (2025-05-22)

* Cache fixes
* Implemented routing

## 1.0.2 (2025-04-01)

* Added 'Protocols' and 'Collections' links to initial page

## 1.0.1 (2025-04-01)

* Fixed error when receiving empty result

## 1.0.0 (2025-04-01)

Integration with [CDD Vault](https://www.collaborativedrug.com/cdd-informatics-platform) registration system.
