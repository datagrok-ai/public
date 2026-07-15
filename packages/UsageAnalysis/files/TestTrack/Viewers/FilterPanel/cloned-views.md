---
feature: filters
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: []
realized_as:
  - cloned-views-spec.ts
related_bugs: []
---

1. Open spgi-100
2. Open the **Filter Panel**
3. On the Filter Panel, navigate to Competition assay filter
4. For Competition assay filter, click the filter menu icon > Missing values > Filter out missing values.
5. On the Filter Panel, set Stereo Category filter to S_ACHIR
6. On the Filter Panel, sketch the c1ccccc1 structure
7. On the Top Menu, go to **View > Layout > Clone View** - a cloned view (spgi-100 copy) opens
8. Check that the FIlter Panel is open and the filtered state on spgi-100 copy matches spgi-100 (Competition assay missing values filtered out, Stereo Category set to S_ACHIR, Structure filter shows c1ccccc1)
9. On the Filter Panel, turn all the filters off (disable the Turn filters on/off checkbox) - check that all rows are visible
10. Turn on filters again (enable the Turn filters on/off checkbox) - check that rows are filtered 
11. On the Filter Panel, locate the Structure filter and click CLEAR, then set the **Structure** filter to C1CCCCC1 - check that filtered state changed
12. On the Filter Panel remove the **Structure** filter (click X icon on its header) - the filtered state should not change (because on the original view the filter should remain present on the Filter Panel)
13. Save the layout
14. Close the Filter Panel
15. Apply the saved layout - the Filter Panel should be open without the Structure filter and the filter state should be the same as before closing the Filter Panel

---
{
"order": 3,
"datasets": ["System:AppData/Chem/tests/spgi-100.csv"]
}
