import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export type FileLoadStategies = {
  demo: (filepath: string) => Promise<DG.DataFrame>,
  dataFiles: (filepath: string) => Promise<DG.DataFrame>,
  eval: (filepath: string) => Promise<DG.DataFrame>,
}

export async function openTableFromDemo(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.getDemoTable(filepath);
  return df;
}
export async function openTablefromDataFiles(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.files.openTable(`System.DemoFiles/${filepath}`);
  return df;
}

export async function openTableWithEval(filepath: string): Promise<DG.DataFrame> {
  const [df] = await (grok.functions.eval(`OpenServerFile("System:DemoFiles/${filepath}")`));
  grok.shell.addTableView(df);
  return df;
}

export const fileLoadStategies: FileLoadStategies = {
  demo: openTableFromDemo,
  dataFiles: openTablefromDataFiles,
  eval: openTableWithEval,
};

export async function addTablesFromPacakgeCsv(_package: DG.Package): Promise<void> {
  // Recursively list package files
  const files = await _package.files.list('', true);

  // Filter files by extension
  const csvFiles = files.filter((f) => f.extension === 'csv');

  // Load every table and add a view for it
  for (const file of csvFiles) {
    // Alternative ways to read a table are:
    // const df = await grok.data.loadTable(`${_package.webRoot}${file.name}`);

    // const df = await _package.files.readCsv(file.name);
    const df = await grok.data.files.openTable(`System:AppData/${_package.name}/${file.fileName}`);
    grok.shell.addTableView(df);
  }
}
