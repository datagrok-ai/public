/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './dataframe/df';
import {tests} from "./test";

export let _package = new DG.Package();

//name: test
export async function test() {

  let results = tests.map(async (t) => {
    let r: {category?: String, name?: String, success: boolean, result: String};
    try {
      r = {success: true, result: await t.test()};
    } catch (x: any) {
      r = {success: false, result: x.toString()};
    }
    r.category = t.category;
    r.name = t.name;
    return r;
  });

  let data = await Promise.all(results);
  grok.shell.addTableView(DG.DataFrame.fromObjects(data)!);
}
