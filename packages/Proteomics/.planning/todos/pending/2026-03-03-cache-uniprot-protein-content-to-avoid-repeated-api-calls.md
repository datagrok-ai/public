---
created: 2026-03-03T17:14:00.000Z
title: Cache UniProt protein content to avoid repeated API calls
area: analysis
files:
  - packages/Proteomics/src/analysis/ (future UniProt panel code)
---

## Problem

When a user clicks on different proteins and then revisits one, the UniProt API is called again for data that was already fetched. This adds latency and unnecessary network traffic, especially when browsing through a protein list.

## Solution

Add a client-side cache (e.g. `Map<string, UniProtEntry>`) keyed by protein/accession ID. On context panel request, check cache first before calling the UniProt REST API. Consider a reasonable cache size limit (LRU) or session-scoped lifetime so it doesn't grow unbounded.
