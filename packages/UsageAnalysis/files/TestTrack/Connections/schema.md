---
feature: connections
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: []
realizes: [views.connections]
realized_as:
  - 07-schema.test.ts
related_bugs: []
---

1. Go to **Browse > Databases**.
2. Expand the **Postgres**.
3. Right-click the `Northwind` connection and select **Browse** from the context menu.
4. Check that you can interact with the structures presented in schema view as with DB Tables (from its context menu). 

---
{
"order": 6
}
