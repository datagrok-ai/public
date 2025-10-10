### Create metadata schema & entity type (by Admin)

1. Log in with admin permissions.
2. Navigate to Browse → Platform → Sticky Meta → **Types**
3. Click "New Entity Type..." button.
4. Enter name (e.g. **TestEntity1**), set matching expression (e.g. semtype=molecule). Save.
5. Navigate to Browse → Platform → Sticky Meta → **Schemas**
6. Click "NEW SCHEMA" button.
7. Enter schema name (e.g. **TestSchema1**), associate with created metadata object type (e.g. "TestEntity1").
8. Add properties: rating (int), notes (string), verified (bool), review_date (datetime). Save schema.

- Expected Result: 
  - Schema and type are created and listed in respective sections.
  - Schema appears in Schemas list; Type appears in Types list; no errors.


---
{
  "order": 1
}
