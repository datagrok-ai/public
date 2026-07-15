---
feature: queries
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: []
realized_as:
  - postgres-query-lifecycle.test.ts
related_bugs: []
---

1. Go to **Browse** > **Databases** > **Postgres** > **NorthwindTest** and find query from the previous steps with name `new_test_query`
2. Delete `new_test_query`:
    * Right-click the connection and select **Delete** from the context menu
    * In the confirmation dialog, click DELETE
3. Reftesh **Browse** view - verify that query has been deleted and is no longer present
---
{
  "order": 4
}