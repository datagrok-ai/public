1. Go to Data > Databases > Postgres 
1. Right-click the `test_postgres` connection
1. Select **Configure Identifiers...** from the context menu
1. Set Schema to `public` and click OK
1. Add a new identifier:
   - Semantic Type: `CUSTOMER_ID`
   - Table: `customers`
   - Column: `customerid`
   -  MatchRegexp: `[A-Z]{5}`
1. SAVE 
1. Reload the page
1. Go to Data > Databases > test_postgres > Schemas > public
1. Open the `customers` table - 
**Values in the customerid column should be highlighted in blue**
1. Click the column's header and check the semantic type (Context Panel > Details)
1. Go to the `test_postgres` connection and remove the identifiers configuration (context menu > Remove)
1. Reload the page 
1. Verify that customerid column values are not highlighted any more
---

{
"order": 1
}
