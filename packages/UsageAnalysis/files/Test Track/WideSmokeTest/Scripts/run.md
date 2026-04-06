1. Go to **Browse > Platform > Functions > Scripts**.
2. Find the `testRscript` script from previous test and right-click it.
3. Select **Run** from the context menu.
4. In the parameters dialog,choose the sample dataset (cars), and click **OK**.
4. Rerun the script, choose any dataset from the local machine, and click **OK**.
4. Rerun script, choose any dataset from Datagrok Files using folder icon, and click **OK**.
4. Rerun script, choose any query using the datasource icon, and click **OK**.
5. Open Datagrok console (Alt + C).
6. Enter the following expression into the console: *{your namespace}:testRscript("cars")*. Where {your namespace}
   equals your platform login.
7. Press enter to execute the command typed in the console.
8. You should get green output in the console with the script output

---
{
"order": 3
}
