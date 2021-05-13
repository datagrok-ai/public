/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import {LaunchBrowser} from "./main.js"

export let _package = new DG.Package();

//name: Repertoire Browser
//tags: app
export async function RepertoireBrowserApp() {

    // let loaded;
    // if (loaded == undefined) {
    //     let tname = grok.shell.tableNames;
    //     let view = null;
    //     if (tname === null || tname.length === 0) {
    //         let df = (await grok.functions.eval('OpenServerFile("Dskatov:RepertoireBrowser/RepertoireBrowserSample.csv")'))[0];
    //         view = grok.shell.addTableView(df);
    //     } else {
    //         view = grok.shell.getTableView(tname[0]);
    //     }
    //     grok.shell.v = view;
    //     await launchBrowser(view);
    //     loaded = true;
    // }

    let tname = grok.shell.tableNames;
    let view = grok.shell.getTableView(tname[0]);
    grok.shell.v = view;

    let app = new LaunchBrowser();
    await app.init(view);
}
