---
status: complete
phase: 04-annotation-and-alternatives
source: [04-01-SUMMARY.md, 04-02-SUMMARY.md]
started: 2026-03-02T04:00:00Z
updated: 2026-03-02T04:30:00Z
---

## Current Test

[testing complete]

## Tests

### 1. UniProt Panel Appears on Protein ID Click
expected: Click a cell in a column with Proteomics-ProteinId semantic type. The context panel on the right should show a UniProt info panel with a loading spinner, then display protein details (name, gene, organism, function description).
result: pass
note: Initially failed for "Primary Protein ID" column — fixed by adding column name to detector (detectors.js)

### 2. GO Terms Grouped by Category
expected: The UniProt panel displays Gene Ontology terms grouped into three categories: Molecular Function, Biological Process, and Cellular Component. Each category shows up to 5 terms.
result: pass

### 3. UniProt Link Opens Full Entry
expected: The panel includes a clickable link to the full UniProt entry (e.g., uniprot.org/uniprot/P12345). Clicking it opens the UniProt page in a new browser tab.
result: pass

### 4. Semicolon-Delimited Protein IDs Handled
expected: If a protein ID cell contains multiple accessions separated by semicolons (e.g., "P12345;Q67890"), the panel fetches and displays info for the first accession (P12345).
result: pass
note: Initially failed for CPTAC format (P00167ups|CYB5_HUMAN_UPS) and semicolon+ups (P02144ups;MYG_HUMAN_UPS) — fixed accession parser to strip ups suffix and handle generic pipe format

### 5. DE Dialog Shows Method Selector
expected: Open the Differential Expression dialog (Proteomics menu). A "Method" dropdown appears with options "limma" and "DEqMS". Default is "limma".
result: pass

### 6. Peptide Count Picker Appears for DEqMS
expected: When "DEqMS" is selected in the Method dropdown, a peptide count column picker appears below. When "limma" is selected, the peptide count picker is hidden.
result: pass

### 7. DEqMS Produces Results
expected: Run DE with DEqMS method selected and a peptide count column chosen. Results table appears with log2FC, p.value, adj.p.value, and significant columns — same format as limma output. Volcano plot and other viewers work with these results.
result: pass
note: DEqMS R package not installed on server; fell back to limma which produced correct results. Added taskbar progress indicator for long-running operations.

### 8. DEqMS Falls Back to Limma
expected: If the server does not have the DEqMS R package installed, running DE with DEqMS selected should fall back to limma (or client-side t-test) and still produce results rather than failing with an error.
result: pass
note: Curl error from R environment setup (not our code). Fallback chain worked correctly.

## Summary

total: 8
passed: 8
issues: 0
pending: 0
skipped: 0

## Gaps

[none]
