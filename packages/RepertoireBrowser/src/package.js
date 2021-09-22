/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { LaunchBrowser } from './main.js';
import { TreeBrowser } from './tree.js';

export let _package = new DG.Package();

//name: Repertoire Browser
//tags: app
export async function RepertoireBrowserApp() {
  let tnames = grok.shell.tableNames;
  let view = null;

  if (tnames === null || tnames.length === 0) {
    try {
      //let df = (await grok.functions.eval('OpenServerFile("Dskatov:RepertoireBrowser/RepertoireBrowserSample.csv")'))[0];
      let df = DG.DataFrame.create(5);
      view = grok.shell.addTableView(df);
    } catch (e) {
      grok.shell.warning('File "Dskatov:RepertoireBrowser/RepertoireBrowserSample.csv" not found.');
      console.log(e);
      return;
    }
  } else {
    view = grok.shell.getTableView(tnames[0]);
  }
  grok.shell.v = view;

  let app = new LaunchBrowser();
  await app.init(view);

  // let df = grok.shell.table('First500seqs');
  // if (df) {
  //   let treeBrowser = new TreeBrowser();
  //   await treeBrowser.init(df);
  // }
}
