---
feature: scripts
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: []
realized_as:
  - scripts-edit-debugged.test.ts
related_bugs: []
---

1. Go to **Browse > Platform > Functions > Scripts**
2. Find the `testRscript` script from Create test case and double-click it
3. Add the expression `newParam="test"` to the script body
4. On the toolbar, click the**Save** button
5. Click on the **x** icon to close script view
6. Double-click `testRscript` script again to open it and check that `newParam="test"` is present in the script

---
{
"order": 2
}
