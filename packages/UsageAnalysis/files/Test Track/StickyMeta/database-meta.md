### Database meta

- Precondition: Context panel is visible
- Known issues:
    - GROK-19427: Database meta: Columns: error on saving: 'Value test doesn't match expected type string_list'
    - GROK-19429: Database meta: Row count deleting is not working

#### Database metadata display
- Open Browse > Databases > Postgres > **CHEMBL**.
- Observe the Context Panel. 
    - Expected: **Database meta** section is present.
- **Fill** database metadata:
    - Comment: test
    - LL M Comment: test
- Click "Save" button.    
- Reload platform. Open Browse > Databases > Postgres > CHEMBL again.
    - Expected: Previously entered metadata is saved and displayed correctly.
- **Clear** all previously added data. Save. Reload page. 
    - Expected: Comment and LLM Comment fields are empty. Metadata cleanup does not affect other tables. No residual metadata appears on other tables.

#### Table metadata display
- Open Browse > Databases > Postgres > **NorthwindTest > Schemas > Public**
- Observe the Context Panel. 
    - Expected: **Database meta** section is present.
- **Fill** database metadata:
    - Comment: test@#$!
    - LL M Comment: test@#$!
- Press "Save" button.    
- Reload platform. Open Browse > Databases > Postgres > CHEMBL.
    - Expected: Previously entered metadata is saved and displayed correctly.
- **Clear** all previously added data. Save. Reload page. 
    - Expected: Comment and LLM Comment fields are empty. Metadata cleanup does not affect other tables. No residual metadata appears on other tables.

#### Column metadata display
- Open Browse > Databases > Postgres > **NorthwindTest > Schemas > Public > categories > categoryid** column.
- Observe the Context Panel. 
    - Expected: **Database meta** section is present.
- **Fill** all the database metadata fields and save.
- Reload platform and open same column.
    - Expected: All previously entered metadata fields are saved and displayed correctly.
- **Clear** all previously added metadata. Save. Reload page. 
    - Expected: All metadata fields are empty. Is unique checkbox is unchecked Clearing metadata does not affect other columns or tables.  

#### Notes
- Pay special attention to known issues when validating save behavior (GROK-19427 and GROK-19429).
- No metadata leakage should occur between:
    - databases,
    - tables,
    - columns.

---
{
  "order": 4
}

