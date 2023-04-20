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

const validFileNames = ['1bdq.pdb', '1bdq.sdf', 'dc.mol2', '4tkx.mmcif', ''];

category("MolstarPreview", () => {
  test("openpdb", async () => {
    let noException = true;
    const pdbFolder: string = `System:AppData/${_packageName}/samples`;
    const pdbFile = (
      await grok.dapi.files.list(pdbFolder, false, 'mol1k.sdf')
    )[0];

    try {
      const view = molecule3dNglView1(pdbFile);
      grok.shell.newView("Molstar Preview", [view]);
    } catch (e) {
      noException = false;
    }
    expect(noException, true);
    await new Promise(r => setTimeout(r, 4000));
  });
});
