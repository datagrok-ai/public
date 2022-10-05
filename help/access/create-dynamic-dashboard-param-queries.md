# Create dynamic dashboards for query results

Follow the steps below to convert a parameterized query output into a dynamic dashboard.

Prerequisites:

1. Connect a database using supported connectors. For instructions,
   see [Add new connecition](databases.md/#add-new-connection).
2. Create a parameterized query. For instructions, see [Parameterized queries](databases.md/#parameterized-queries).

Create a dynamic dashboard:

1. Go to **Data** > **Databases**.
2. In the **Database Explorer**, right-click the desired query and select **Run**. A parameter dialog opens.
3. Fill in the query parameters and click **OK**. A dataframe opens.
   > Note: A query output has a unique URL that encodes the parameters entered. Use the **Toolbox** on the left to
   change the parameters and note the URL change. You can share the query results using this link.
4. Optional. Add viewers, filter, or apply layouts to customize the dashboard.
   > Tip: You can save a _layout_ as a template for future use. To do so, on the **Toolbox** on the left, expand the **
   Layout** info panel and click **Save**. The layout appears in the gallery. To learn more about _layouts_,
   see [Layouts](../visualize/view-layout.md).
5. Once you have the dashboard set up the way you want it, save it by uploading athe project to the server:
   1. On the **Sidebar** click **Projects** and then click the **Upload** button on the **Toolbox**. A dialog opens.
   2. In the dialog, enter the name for the project and, optionally, the project description in the fields provided.
   3. Select how to store data. You have two options: (1) save the data as a static snapshot, and (2) store the data as a
      generation script. When you select the second option, the query is executed every time you open the project. To
      select this option, toggle the **Data sync** control.
      > Note: To learn more about dynamic data updates in projects,
      > see [Dynamic data](../datagrok/project.md/#dynamic-data).
6. Click **OK** to upload the project.

Once the projects has been uploaded to the server, you can open it or share it with others. Go to **Data** > **
Projects** and locate your project. Right-clicking the project opens its context menu, double-clicking the project opens
it in Datagrok.

![Dynamic dashboards](dynamic-dashboards.gif)
