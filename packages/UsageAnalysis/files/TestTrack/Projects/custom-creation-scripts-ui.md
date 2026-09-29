---
feature: projects
realizes_atlas: [projects.cp.upload-save-reopen-golden]
realizes: [views.projects]
priority: p0
target_layer: playwright
coverage_type: regression
related_bugs: []
---

# Custom creation script — the project picks up the newest file on reopen

A project's table comes from your own script. The script reads a folder
in **My files** and returns the CSV file with the highest number at the
end of its name. The project is saved with **Data sync**, a file with a
higher number is added to the folder, and the reopened project must
show the new file's rows, not the rows it was saved with.

## Setup

1. Log in as the test user.
2. Names in this test: folder `ccsTest`, script `ccsScript`, project
   `ccsProject`.
3. **Create the folder.**
   - Go to **Browse > Files**.
   - Right-click **My files** and choose **Create folder...**.
   - Replace *Folder name* with `ccsTest`.
   - Click **OK**.
4. **Put the first file into it.**
   - Go to **Browse > Files > Demo**.
   - Drag `cars.csv` onto the `ccsTest` folder under **My files** in
     the Browse tree.
   - In **My files > ccsTest**, right-click `cars.csv` and choose
     **Rename...**.
   - Change the name to `data1.csv`.
   - Click **OK**.
   - **Verify:** `ccsTest` holds only `data1.csv`.
5. **Create the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Click **NEW** and choose **JavaScript Script...**.
   - Replace the template with the text below.
   - Click **SAVE** in the editor.
   - **Verify:** the balloon *Script saved.* appears.

   ```js
   //name: ccsScript
   //language: javascript
   //output: dataframe result
   const folder = grok.shell.user.project.name + ':Home/ccsTest/';
   const csvFiles = await grok.dapi.files.list(folder, false, 'csv');
   if (csvFiles.length === 0)
     throw new Error('No CSV files found in ' + folder);
   const suffix = (name) => { const m = name.match(/(\d+)(?=\.csv$)/); return m ? parseInt(m[1], 10) : -1; };
   csvFiles.sort((a, b) => suffix(a.fileName) - suffix(b.fileName));
   result = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(csvFiles[csvFiles.length - 1].fullPath));
   ```

## Scenario

1. **Run the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Click the refresh icon.
   - Right-click `ccsScript` and choose **Run...**.
   - **Verify:** the `result` view opens with **Rows: 30** in the status
     bar (the rows of `data1.csv`).

2. **Add a viewer and save the project.**
   - In **Toolbox > Viewers**, click **Bar chart**.
   - Click **SAVE** on the ribbon.
   - **Verify:** **Data sync** is ON for `result`.
   - Expand **CREATION SCRIPT**.
   - **Verify:** it shows a call of `ccsScript`.
   - Enter `ccsProject` as the name.
   - Click **OK**.
   - **Verify:** the balloon *Project "ccsProject" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.

3. **Close everything.**
   - Right-click the left sidebar and select **Close All**.
   - **Verify:** no table view is open.

4. **Add a file with a higher number.**
   - Go to **Browse > Files > Demo**.
   - Drag `iris.csv` onto the `ccsTest` folder under **My files**.
   - In **My files > ccsTest**, right-click `iris.csv` and choose
     **Rename...**.
   - Change the name to `data2.csv`.
   - Click **OK**.
   - **Verify:** `ccsTest` holds `data1.csv` and `data2.csv`.

5. **Reopen the project.**
   - Go to **Browse > Dashboards**.
   - Type `ccsProject` into the search box.
   - Click the refresh icon.
   - Double-click the `ccsProject` tile.
   - **Verify:** the `result` view opens with **Rows: 150** (the rows of
     `data2.csv`, not the 30 it was saved with).
   - **Verify:** the view has the bar chart.
   - **Verify:** no error balloon and no **Data loading error** dialog
     appear.

## Cleanup

- Right-click the left sidebar and select **Close All**.
- In **Browse > Dashboards**, right-click the `ccsProject` tile and
  choose **Delete Project**. Click **DELETE** and wait until the dialog
  closes.
- In **Scripts**, right-click `ccsScript` and choose **Delete**. Click
  **YES**.
- In **Browse > Files > My files**, right-click `ccsTest` and choose
  **Delete...**. Confirm.

## Expected results

- A table produced by your own script is saved with the script as its
  creation script.
- On reopen the script runs again and the table shows what the script
  returns now: the newest file's rows.

## Automation notes

- Copying a demo file into **My files** by dragging it onto the folder
  in the Browse tree is read from `files_view.dart` (a drop onto a
  folder of another connection copies the file); it was not run by hand
  while writing this file. The two files are setup, not the subject, so
  an automated run may write them through the JS API.
- The confirmation text of **Delete...** for a folder was not checked
  in source.

---
{
  "order": 2,
  "datasets": []
}
