1. Reload the tree in Browser
2. Right-click the connection from the previous step (Adding).
3. Select **Edit** from context menu.
4. Change name to `new_test_postgres` and click OK.
5. Change login/password parameters of the `new_test_postgres` connection with arbitrary data and save the changes.
6. Test the connection - it should return an error:

   ```
   "new_test_postgres": failed to connect:
   ERROR:
   com.zaxxer.hikari.pool.HikariPool$PoolInitializationException: Failed to initialize pool: FATAL: password authentication failed for user "new_test_postgres"
   ```

7. Set the right login/password parameters and test the connection - it should be fine.

* Repeat scenarios 2.1-2.2 for connections to different DB sources (Oracle, MariaDB, MySQL, MS SQL)
---
{
  "order": 2
}
