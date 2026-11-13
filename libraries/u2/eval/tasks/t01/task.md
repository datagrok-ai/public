# Submit a campaign

When a Hit Design campaign is ready, the chemist submits it from a **Submit campaign** dialog that sums the campaign
up first.

- Rows **Name**, **File Path** and **Status** show the campaign's values, read-only: the user can select and copy
  them.
- The Name row has a pencil icon, **Edit name**, for renaming the campaign in place. It puts the cursor in Name;
  Enter keeps the new name (`actions.rename(newName)`), Escape brings the old one back. Either way Name is
  read-only again afterwards.
- The File Path row has a **Download** icon that calls `actions.download()`.
- **SUBMIT** calls `actions.submit()`. CANCEL closes the dialog.

The plugin calls this, and it shows the dialog and returns it:

```ts
export function submitDialog(campaign: {name: string, filePath: string, status: string},
  actions: {rename(name: string): unknown, download(): void, submit(): unknown}): Dialog
```
