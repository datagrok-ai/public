---
title: "Building light-weighted hit triage system"
---

With Datagrok, creating a light-weighted hit triage system using your database as source is easy and efficient. By applying main Datagrok features, you can effectively screen molecules.

- Query database

  ```js
  let data = await grok.functions.call('Chembl:_compoundActivityDetailsForTarget', {target: "CHEMBL1827"})
  ```

Here we used an existing data query. However, if you want to learn more about query creation, go to this [link](https://datagrok.ai/help/develop/how-to/access-data#creating-queries).

- Add required columns to the received results and add tableView

  ```js
  data.columns.addNewString('Comments')
  ```

- Set column edit privileges if required

    ```js
    data.col('Comments').setTag('editableBy', 'askalkin, aparamonov')
    ```

For more information on setting privileges, visit this [link](https://datagrok.ai/help/visualize/viewers/grid#column-edit-permissions)

- Add data as table view

    ```js
    await grok.data.detectSemanticTypes(data)
    let tableView = grok.shell.addTableView(data)
    ```

- Apply layout if required

    ```js
    let layouts = await grok.dapi.layouts.getApplicable(data)
    console.log(layouts)
    tableView.loadLayout(layouts[0]);
    ```

For more information on manipulating layouts, check out this [link](https://datagrok.ai/help/develop/how-to/layouts)
