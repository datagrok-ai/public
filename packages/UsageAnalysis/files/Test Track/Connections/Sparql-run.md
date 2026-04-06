# Sparql — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Press three dots at end of Databases list | PASS | PASSED | The full provider list is visible after scrolling — Sparql appears after Snowflake; the "..." expands hidden providers; Sparql node found in tree |
| 2 | Right-click Sparql and select Add new connection | PASS | PASSED | Right-clicked Sparql tree node; menu: "Browse queries", "Browse connections", "New connection..."; clicked "New connection..." |
| 3 | Enter `test_sparql` to the Name field | PASS | PASSED | "Add new connection" dialog for Sparql type appeared with Name, Endpoint, Requires Server, Prefixes fields; Name set to test_sparql |
| 4 | Fill Endpoint, Requires Server, Prefixes | PASS | PASSED | Endpoint=http://data.ontotext.com/repositories/data-last, Requires Server=checked, Prefixes=empty |
| 5 | Click the Test button | FAIL | FAILED | Test returned: `"test_sparql": failed to connect: SocketException: Failed host lookup: 'data.ontotext.com' (OS Error: Name or service not known, errno = -2)` — endpoint hostname not resolvable from server |
| 6 | Click OK | PASS | PASSED | OK clicked; dialog closed; test_sparql saved and appears under Sparql node in tree |
| 7 | Delete test_sparql connection | PASS | PASSED | Right-clicked test_sparql → Delete... → DELETE confirmed; node removed from tree |

## Summary

6 of 7 steps passed (1 failed). All UI steps worked correctly. The SPARQL connection was created and deleted successfully. Step 5 failed because the test endpoint `data.ontotext.com` is not resolvable from the Datagrok server — the endpoint may be discontinued or blocked.

## Retrospective

### What worked well
- Sparql provider is accessible via the expanded "..." section of the Databases tree
- "New connection..." context menu option creates a SPARQL-specific dialog
- SPARQL dialog has the correct field layout: Name, Endpoint, Requires Server checkbox, Prefixes textarea
- Connection saves and deletes correctly regardless of test result

### What did not work
- `data.ontotext.com` endpoint not reachable: `SocketException: Failed host lookup` — endpoint appears to be down or blocked from the server

### Suggestions for the platform
- Provide a reliable public SPARQL test endpoint in the scenario documentation
- Show connection type label in the "Add new connection" dialog title (e.g., "Add new Sparql connection")

### Suggestions for the scenario
- Step 5: Update with a currently working SPARQL endpoint for testing (data.ontotext.com appears to be unavailable)
- Consider using `https://query.wikidata.org/sparql` or another reliable public endpoint
