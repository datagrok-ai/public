import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { MolecularLiabilityBrowser } from './molecular-liability-browser';
//import { TreeBrowser } from './tree.js';

export let _package = new DG.Package();

function getPathSegments(path: string) {
  let parser = document.createElement('a');
  parser.href = path;
  let pathSegments = parser.pathname.split('/');
  if (pathSegments.length > 4)
    return pathSegments[4];
  else
    return null;
}

//name: Repertoire Browser
//tags: app
export async function RepertoireBrowserApp() {
  grok.shell.windows.showToolbox = false;
  let vid = getPathSegments(<string><unknown>window.location);

  let pi = DG.TaskBarProgressIndicator.create('Opening Molecular Liability Browser');
  let app = new MolecularLiabilityBrowser();
  await app.init(vid);
  pi.close();
  //let tnames = grok.shell.tableNames;
  //let view = null;

  // if (tnames === null || tnames.length === 0) {
  //   try {
  //     //let df = (await grok.functions.eval('OpenServerFile("Dskatov:RepertoireBrowser/RepertoireBrowserSample.csv")'))[0];
  //     const path = _package.webRoot + 'tableFiles/mlb.csv';
  //     const df = (await grok.data.loadTable(path));
      

  //     view = grok.shell.addTableView(df);
  //   } catch (e) {
  //     grok.shell.warning('File "Dskatov:RepertoireBrowser/RepertoireBrowserSample.csv" not found.');
  //     console.log(e);
  //     return;
  //   }
  // } else {
  //   view = grok.shell.getTableView(tnames[0]);
  // }

  // let df = grok.shell.table('First500seqs');
  // if (df) {
  //   let treeBrowser = new TreeBrowser();
  //   await treeBrowser.init(df);
  // }
}
