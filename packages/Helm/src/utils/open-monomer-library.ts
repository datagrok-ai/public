import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export async function openMonomerLibrary(fi: DG.FileInfo): Promise<void> {
  const obj = JSON.parse(await fi.readAsString());
  const df = DG.DataFrame.fromObjects(obj);
  if (df)
    grok.shell.addTableView(df);
  else
    grok.shell.warning('Failed to open monomer library');
}
