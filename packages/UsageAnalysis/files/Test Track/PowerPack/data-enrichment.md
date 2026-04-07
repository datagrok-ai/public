## 1. Create, edit, and remove enrichment

1. Go to Databases > Postgres > Northwind
1. Create a SQL and a visual query for the 'Orders' table
1. Open the Orders table and run created queries
2. In the opened table, click the 'customerid' column
1. Go to the **Context Panel > Northwind** (or other connection name you use) **> Enrich..**
1. Click **+ Add Enrichment..**
3. Click **Add a table to join > public > Customers** and select some columns 
1. Verify the editor shows joins correctly
1. Enter an enrichment name and click **SAVE** 
4. Click the created enrichment - **the selected columns should be added to the grid**
1. Edit the enrichment by changing selected columns or join settings

## 2. Multiple enrichments per column and multiple columns

1. On the same column, create a second enrichment with different columns
2. Create an enrichment for a different column in the same table
3. Apply all enrichments
1. Remove any enrichment

## 3. Persistence across projects, layouts, and reuse in other tables

1. Go to the query results - verify that created enrichments are available there as well
1. Create new enrichments for any columns
1. Apply enrichments and save the project/layout
1. Delete the added columns and apply the layout - verify that enriched columns are present again
2. Close and reopen the project
3. Open another table (e.g. customercustomerdemo) or query result with the same ('customerid') column
4. Verify that previously created enrichments are available for reuse

## 4. Visibility for different users

1. Log in with another user
1. Go to Databases > Postgres > Northwind (or another connection you used)
1. Open the Orders table
1. Click the 'customerid' column
3. Verify that previously created enrichments are available for reuse

---
{
  "order": 4
}