### External Provider

#### Navigate to PostgreSQL Connection:
1. Browse > Databases > Postgress.
2. Click on **Add connection...**. 
3. Fill out the "Add new connection" form with the following details:
* Data Source: Postgres
* Name: PostgreSQLDBTests2
* Server: db.datagrok.ai
* Port: 54327
* Db: test
* Login: superuser
* Password: *** (obtain from QA or DevOps)

#### Create and Run Queries:
4. For the created connection, add and run four queries sequentially. 

* **TestCreateTable**. Save and run the query:

CREATE TABLE tmp_table_test (id bigint, name varchar);

* **TestInsertData**. Save and run the query:

INSERT INTO tmp_table_test VALUES (1, 'test');

* **TestUpdateData**. Save and run the query:

UPDATE tmp_table_test SET name = 'test' WHERE id = 1;

* **TestDropTable**. Save and run the query:

DROP TABLE tmp_table_test;

5. Delete PostgreSQLDBTests2 connection.

#### Expectes results: 
* The connection should be created and deleted successfully.
* All queries should run successfullys.
* No errors should occur during the test.

---
{
"order": 8
}
