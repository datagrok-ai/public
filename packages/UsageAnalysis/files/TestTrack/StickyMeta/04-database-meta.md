---
feature: stickymeta
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: []
realized_as:
  - 04-database-meta.test.ts
related_bugs: [GROK-19427, GROK-19429]
---

# 04 — Database meta on a DB table and column

Verifies the **Database meta** context-panel section for a database **table** and **column** —
filling, saving, persistence across reload, clearing, and no metadata leakage between objects.

**Preconditions**

- A reachable Postgres connection with an introspectable schema. Recon used **`NorthwindTest`**
  (`Browse > Databases > Postgres > NorthwindTest > Schemas > public`), table **`categories`**,
  column **`categoryid`**.
- Context Panel is visible.

**Scope**

- Only **table** and **column** levels are covered. The **connection (DbInfo)** entity does
  **not** expose a Database meta section in the UI on dev, so connection-level metadata is out of
  scope for this scenario.

**Known issues to watch** (validate save behavior against these)

- **GROK-19427** — Database meta > Columns: error on saving `Value test doesn't match expected
  type string_list`.
- **GROK-19429** — Database meta: Row count deleting is not working.

> Navigation: expand `Databases > Postgres > NorthwindTest`, then the **Schemas** node →
> `public` → `categories` → `categoryid`. Selecting a table/column node shows its entity in the
> Context Panel with a **Database meta** accordion section.

---

## Test 4.1 — Table metadata display, save, persist, clear

**Steps**

1. Open `Browse > Databases > Postgres > NorthwindTest > Schemas > public` and select the
   **`categories`** table.
2. In the Context Panel, expand **Database meta**. Confirm the section is present with fields
   **Domains**, **Row Count**, **Comment**, **LLM Comment** and a **Save** button.
3. Fill `Comment` = "test" and `LLM Comment` = "test". Click **Save**.
4. **Reload** the platform, navigate back to `categories`, expand Database meta.
5. **Clear** Comment and LLM Comment. Click **Save**. **Reload** and re-open.

**Expected**

- The Database meta section is present on the table.
- After step 3 + reload: `Comment` = "test" and `LLM Comment` = "test" persist.
- After step 5 + reload: both fields are empty.
- No metadata appears on any other table (no leakage between tables).

---

## Test 4.2 — Column metadata display, save, persist, clear

**Steps**

1. Expand `categories` and select the **`categoryid`** column.
2. In the Context Panel, expand **Database meta**. Confirm the section is present with fields
   including **Is Unique**, **Min**, **Max**, **Values**, **Sample Values**, **Unique Count**,
   **Quality**, **Comment**, **LLM Comment** and a **Save** button.
3. Fill `Comment` and `LLM Comment` (and optionally toggle **Is Unique**). Click **Save**.
4. **Reload** the platform, re-open the same column, expand Database meta.
5. **Clear** the metadata. Click **Save**. **Reload** and re-open.

**Expected**

- The Database meta section is present on the column.
- After step 3 + reload: the saved values persist (watch GROK-19427 on save).
- After step 5 + reload: the fields are empty; **Is Unique** is unchecked.
- Clearing column metadata does not affect any other column or table (no leakage).

---

## Cleanup

- Ensure any Comment / LLM Comment / Is Unique values set during the test are cleared and saved
  (so the shared connection is left in its original state).
