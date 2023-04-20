import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import {
  category /*, expect*/,
  test,
  expect
} from "@datagrok-libraries/utils/src/test";
import { _package, molecule3dNglView1 } from "../package";
import { _packageName } from "./utils";

const validFileNames = ['1bdq.pdb', '1bdq.sdf', 'dc.mol2',
  '4tkx.mmcif', 'caffeine.xyz', 'grofile.gro', 'pdbqt.pdbqt'];

category("MolstarPreview", () => {

  validFileNames.forEach(fn =>{
    test(`open${fn.substring(fn.indexOf('.'),fn.length)}`, async () => {
      let noException = true;
      const folderName: string = `System:AppData/${_packageName}/samples`;
      const file = (
        await grok.dapi.files.list(folderName, false, fn))[0];

      try {

        const view = molecule3dNglView1(file);
        grok.shell.newView("Molstar Preview", [view]);
      } catch (e) {
        noException = false;
      }
      expect(noException, true);
    });
  });

  test('openCsvFile', async () => {
    let noException = true;
    const folderName: string = `System:AppData/${_packageName}/samples`;
    const file = (
    await grok.dapi.files.list(folderName, false, 'dock.csv'))[0];

    try {
      const view = molecule3dNglView1(file);
      grok.shell.newView("Molstar Preview", [view]);
    } catch (e) {
      noException = false;
    }
    expect(noException, false);
  });

});
