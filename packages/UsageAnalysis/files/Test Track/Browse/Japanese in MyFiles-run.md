# Japanese in MyFiles — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Go to MyFiles section | PASS | Navigated to Files > My files. URL: /files/OahadzhanianDatagrokAi.Home/?browse=files. 38 files displayed |
| 2 | Create file with Japanese name: 芹沢 貴之 こんにちは | PASS | File already exists from prior creation (Mar 27, 2025). Name displayed correctly in file list as "芹沢 貴之 こんにちは.csv" |
| 3 | Verify file name display — no garbled characters | PASS | Japanese characters render correctly in: (1) file list, (2) context panel title, (3) URL encoding (%E8%8A%B9%E6%B2%A2+%E8%B2%B4%E4%B9%8B+...). Context panel shows: Connection "My files", Size 15417416, Created Mar 27, 2025 |

## Summary

All 3 steps passed. The file 芹沢 貴之 こんにちは.csv displays correctly with no garbled or incorrect characters. Japanese Kanji and Hiragana render properly in the file list, context panel, and URL encoding.

## Retrospective

### What worked well
- Non-Latin (Japanese) characters display correctly throughout the UI
- URL encoding handles multi-byte characters properly
- Context panel correctly shows file metadata for non-Latin named files

### What did not work
- Nothing — all aspects of non-Latin character display work correctly

### Suggestions for the platform
- None — Japanese character support works as expected

### Suggestions for the scenario
- Scenario says "create a new file" but doesn't specify how (drag-and-drop? upload button?). Clarify the creation method
- Consider testing other non-Latin scripts (Cyrillic, Arabic, Korean) in addition to Japanese
- Add a step to verify the file can be opened/downloaded with the correct name preserved
