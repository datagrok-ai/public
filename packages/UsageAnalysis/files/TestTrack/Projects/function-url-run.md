---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [views.functions, views.projects]
pyramid_layer: integration
produced_from: manual-draft
related_bugs: [GROK-20511]
---

# Function link — run a function straight from its URL

`/func/<Namespace>.<Name>?<param>=<value>` opens a function's parameter
form. Adding `&run=true` runs it at once and shows only the result. A
copy icon next to **RUN** gives such a link for the current values. A
run link without a required value falls back to the parameter form and
says what is missing.

## Setup

1. Log in as the test user.
2. **Create the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Click **NEW** and choose **JavaScript Script...**.
   - Replace the template with the text below and click **SAVE**.

   ```
   //name: urlRunScript
   //language: javascript
   //input: int rows = 10
   //output: dataframe df
   df = grok.data.demo.demog(rows);
   ```

3. Create the script `urlRunRequired` the same way, with this text:

   ```
   //name: urlRunRequired
   //language: javascript
   //input: int rows
   //input: string label = "x"
   //output: dataframe df
   df = grok.data.demo.demog(rows);
   ```

4. Click `urlRunScript` in **Scripts**. In the **Context Panel**, under
   **Details**, click **Links...** and note the **Grok name**, for
   example `Admin:urlRunScript`. The URL name is the same with a dot:
   `Admin.urlRunScript`. Below, `<server>` is the server address.

## Scenario

1. **The form opens without run.**
   - In the address bar, open `<server>/func/Admin.urlRunScript?rows=25`.
   - **Verify:** the parameter form opens with **Rows** = 25 and a
     **RUN** button; no table is shown yet.

2. **Copy a run link.**
   - Set **Rows** to 40.
   - Hover the copy icon left of **RUN**.
   - **Verify:** the tooltip reads *Copy a link that runs this function
     with the current parameters*.
   - Click the copy icon.
   - Open a new tab, paste, and read the address before pressing Enter.
   - **Verify:** it ends with `/func/Admin.urlRunScript?rows=40&run=true`.

3. **Run from the link.**
   - Press Enter.
   - **Verify:** a table view opens directly, with no parameter form,
     and the status bar shows **Rows: 40**.
   - **Verify:** the address bar still contains `run=true`.
   - **Verify:** **Toolbox > Source** shows **Rows** = 40 and
     **REFRESH**.

4. **Change the value in the result view.**
   - In **Source**, set **Rows** to 15 and click **REFRESH**.
   - **Verify:** the status bar shows **Rows: 15**.
   - Hover the copy icon next to **REFRESH**.
   - **Verify:** the tooltip link ends with `?rows=15&run=true`.

5. **run=false keeps the form.**
   - Open `<server>/func/Admin.urlRunScript?rows=25&run=false`.
   - **Verify:** the parameter form opens; nothing runs.

6. **A required value is missing.**
   - Open `<server>/func/Admin.urlRunRequired?label=y&run=true`.
   - **Verify:** a warning balloon starting with *Unable to run
     "urlRunRequired"* appears and names `rows`.
   - **Verify:** the parameter form opens with **Label** = `y`; nothing
     runs.
   - Enter **Rows** = 5 and click **RUN**.
   - **Verify:** a table with 5 rows opens.

## Cleanup

- Right-click the left sidebar and select **Close All**. In **Scripts**,
  right-click `urlRunScript`, choose **Delete**, click **YES**. Repeat
  for `urlRunRequired`.

## Expected results

- With `run=true` the function runs on open and only the result is
  shown.
- Without it, or with `run=false`, the form opens with the URL values.
- The copy icon gives a run link with the values currently set.
- A run link without a required value falls back to the form and names
  what is missing.

## Automation notes

- Replace `Admin` with the namespace from the script's **Grok name**
  when the test account is not the admin.
- The parameter caption **Rows** (from the input name `rows`) and the
  status-bar reading after **REFRESH** were not checked in source.
- The exact wording after *Unable to run "urlRunRequired":* comes from
  the parameter validation and was not checked in source.
