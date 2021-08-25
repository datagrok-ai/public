/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OpenChemLibSketcher} from "./ocl-sketcher";
import {ProgressIndicator} from "datagrok-api/dg";

export let _package = new DG.Package();

//name: openChemLibSketch
//description: Sketches a molecule
//top-menu: Chem | OpenChemLib Sketch
export function openChemLibSketch() {
  ui.dialog()
    .add(openChemLibSketcher().root)
    .showModal(true);
}

//name: openChemLibSketcher
//tags: moleculeSketcher
//output: widget sketcher
export function openChemLibSketcher() {
  return new OpenChemLibSketcher();
}

//name: openSdf
//description: Opens SDF file
//tags: file-handler
//meta.ext: sdfx
//input: list bytes
//output: list tables
export function openSdf(bytes: Uint8Array) {
  grok.shell.info('opened');
  console.log(bytes.length);
  return [grok.data.demo.demog()];
}