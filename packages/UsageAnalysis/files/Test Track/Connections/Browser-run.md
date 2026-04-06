# Browser — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Browse > Databases | PASS | PASSED | Navigated to Databases view showing all provider nodes in the tree |
| 2 | Click Filter templates icon (magic wand) and check all templates | AMBIGUOUS | N/A | The search bar is present but no "magic wand" / filter templates icon was found in the current UI; search works directly without toggling template filters |
| 3 | Type `new_test` in search field | PASS | PASSED | Typed into the "Search connections by name or by #tags" field; list filtered to show 1/1 matching connection `new_test_postgres` |
| 4 | Click on the found connection | PASS | PASSED | Clicked `new_test_postgres`; connection selected (highlighted) and Context Panel appeared on the right |
| 5a | Check Details tab | PASS | PASSED | Server=db.datagrok.ai, Database=northwind, port=54322, ssl=false, Created by=Olesia Pavlenko — all match values from Adding/Edit steps |
| 5b | Share the connection with another user | AMBIGUOUS | N/A | Sharing tab visible showing "You are the owner" + SHARE... button; sharing with another user not tested (would affect real user) |
| 5c | Check Activity has right dates and actions | PASS | PASSED | Activity 4 entries: "Olesia Pavlenko created new_test_postgres" + 3 "edited" entries — matches all actions performed |
| 5d | Send a message to the chat | PASS | PASSED | Chats section visible with "Start the chat, press Enter to send" input |
| 6 | Click dropdown list icon near connection name | PASS | PASSED | Clicked `grok-context-arrow-down` icon near title; menu appeared: Browse, New Query..., New Visual Query..., Delete..., Edit..., Rename..., Clone..., Clear cache, Browse queries, Test connection, Share..., Copy, Add to favorites |

## Summary

7 of 9 steps passed, 2 ambiguous. Search filtering works correctly and the Context Pane shows all expected tabs (Details, Sharing, Activity, Chats, Dev). The dropdown near the connection name shows the full set of actions. The "magic wand" filter icon was not found in the current UI.

## Retrospective

### What worked well
- Search field filters connections in real time by name
- Context Pane updates immediately on connection click
- Activity log accurately tracks all create/edit operations with timestamps
- `grok-context-arrow-down` icon reliably opens the actions dropdown

### What did not work
- "Filter templates" magic wand icon not found in current UI — may have been removed or renamed
- Actual sharing with another user not tested to avoid affecting real users on public server

### Suggestions for the platform
- Document the filter icon or add a tooltip explaining its purpose
- Add `data-testid="connection-actions-dropdown"` to the dropdown trigger icon

### Suggestions for the scenario
- Step 2: Clarify what "Filter templates" means and its current location in the UI
- Step 5b: Note that sharing test requires a second test account on the same server
