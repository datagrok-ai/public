---
feature: queries
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: []
realizes: [views.queries]
realized_as:
  - postgres-query-lifecycle.test.ts
related_bugs: []
---

1. Right-click `test_query` and select `Edit...` from the context menu
2. Change name to `new_test_query`, and click `SAVE`.
3. Change the query test to: 
```
select * from orders
```
1. Run the query using both:
   * **Menu Ribbon > Play button** — result appears at the bottom of the current view
   * **Toolbox > Actions > Run query…** — result opens in a new view
8. Save the query

---
{
  "order": 2
}