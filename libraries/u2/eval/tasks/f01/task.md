# A campaign's actions in its view

A Hit Design campaign opens in a platform TableView, `tv`. Its actions go on the view's ribbon: the icons **Add
new row** (a plus), **Refresh Campaign** (a sync icon) and **Merge with table**, each with its tooltip, and a
**Submit...** button.

- **Add new row** appends a row to the campaign's table, `tv.dataFrame.rows.addNew(null, true)`. **Refresh
  Campaign** reloads the campaign, `service.refresh()`.
- **Submit...** opens the **Submit campaign** dialog, which sums the campaign up first. Rows **Name**, **File Path**
  and **Status** show the campaign's values, read-only: the user can select and copy them. A pencil on the Name
  row, **Edit name**, renames the campaign in place: it puts the cursor in Name, Enter keeps the new name
  (`service.rename(newName)`), Escape brings the old one back, and either way Name is read-only again
  afterwards. A **Download** icon on the File Path row calls `service.download()`. **SUBMIT** calls
  `service.submit()`.
- **Merge with table** opens the **Merge with table** dialog. Users keep the files they merge on their computer or
  in Datagrok Files, and often paste a file share path such as `System:DemoFiles/demog.csv`; **Source file** takes
  a file from any of the three. Once a file is chosen, `service.readTable(file)` reads it, and **Molecules** asks
  which of its columns holds the structures: it offers the table's molecule columns and starts on the first one it
  offers. **MERGE** calls `service.merge(file, table, moleculeColumn)` and is off until the dialog has what the
  merge needs. A file from Datagrok Files is passed with its share path (`fullPath`), so the merge can be repeated
  later.

The plugin calls this once the view is open:

```ts
type Campaign = {name: string, filePath: string, status: string};

interface CampaignService {
  rename(name: string): Promise<void>;
  download(): void;
  submit(): Promise<void>;
  readTable(file: DG.FileInfo): Promise<DG.DataFrame>;
  merge(file: DG.FileInfo, table: DG.DataFrame, moleculeColumn: string): Promise<void>;
  refresh(): Promise<void>;
}

export function campaignActions(tv: DG.TableView, campaign: Campaign, service: CampaignService): void
```
