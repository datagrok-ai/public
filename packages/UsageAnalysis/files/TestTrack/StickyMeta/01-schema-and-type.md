---
feature: stickymeta
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: []
realizes: []
realized_as:
  - 01-schema-and-type.test.ts
related_bugs: []
---

# 01 — Create / edit / delete entity type and metadata schema

Covers the admin configuration surface of Sticky Meta: defining a metadata **entity type**
(what objects the metadata attaches to) and a metadata **schema** (the typed properties).

**Preconditions**

- Logged in as an **Administrator** (the **NEW ENTITY TYPE...** / **NEW SCHEMA...** buttons and
  delete are admin-only).
- Browse tree is visible.

**Test data**

- Entity type name: `PW_SM_Type_<suffix>`, matching expression `semtype=molecule`.
- Schema name: `PW_SM_Schema_<suffix>`, associated with the entity type above, properties:
  `rating` (int), `notes` (string), `verified` (bool), `review_date` (datetime).

> Navigation note: the Datagrok SPA ignores `goto()` between sibling Sticky-Meta routes
> (`/meta/types` ↔ `/meta/schemas`) — it reloads to the last-visited one. Always switch between
> Types and Schemas by **clicking the tree node**, the way a user does.

---

## Test 1.1 — Create an entity type

**Steps**

1. Go to **Browse > Platform > Sticky Meta > Types**.
2. Click **NEW ENTITY TYPE...**. A dialog "Create a new entity type" opens with fields
   **Name** and **Matching expression**.
3. With both fields empty, observe the **OK** button.
4. Type the name only, observe **OK**.
5. Enter the matching expression `semtype=molecule`. Click **OK**.

**Expected**

- OK is **disabled** while either Name or Matching expression is empty (both are required).
- After OK, the type appears in the Types list (search by name to confirm — the list is paginated).
- No error balloon.

---

## Test 1.2 — Create a schema associated with the type

**Steps**

1. Go to **Browse > Platform > Sticky Meta > Schemas** (click the **Schemas** tree node).
2. Click **NEW SCHEMA...**. A dialog "Create a new schema" opens with **Name**,
   **Associated with:** (a "select entities" link), and a **Properties** table.
3. Enter the schema name.
4. Click **Associated with: > select entities**, check the entity type from Test 1.1 in the
   picker, confirm.
5. Add four properties (use **+** / "Add new property to schema" for each row), setting **Name**
   and **Property Type**:
   - `rating` → `int`
   - `notes` → `string`
   - `verified` → `bool`
   - `review_date` → `datetime`
6. Click **OK**.

**Expected**

- The available property types are `string`, `int`, `bool`, `double`, `datetime`, `string_list`.
- After OK, the schema appears in the Schemas list (search to confirm).
- No error balloon.

---

## Test 1.3 — Edit the schema and verify its content

**Steps**

1. In the Schemas list, search for the schema, right-click its card → **Edit**.
2. The "Edit schema" dialog shows the saved structure.
3. Verify the association and all four properties, then close the dialog.

**Expected**

- **Associated with** shows the entity type from Test 1.1.
- The Properties table lists exactly: `rating`/int, `notes`/string, `verified`/bool,
  `review_date`/datetime, in order.

---

## Test 1.4 — Delete the schema and entity type (cleanup)

**Steps**

1. In the Schemas list, search for the schema, right-click its card → **Delete** → confirm **DELETE**.
2. Switch to **Types** (tree node), search for the entity type, right-click → **Delete** → confirm.

**Expected**

- Both entities disappear from their lists.
- No error balloon.

> **Mandatory cleanup** — Tests 1.1–1.4 form one lifecycle; the schema and type created here
> must be deleted at the end even if an earlier assertion fails (autotest: `finally` block using
> `dapi.stickyMeta.deleteSchema(id)`, see `README.md`).
