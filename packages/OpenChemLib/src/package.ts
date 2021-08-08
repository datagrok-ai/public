/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OpenChemLibSketcher} from "./ocl-sketcher";

export let _package = new DG.Package();

//name: openChemLibSketch
//description: Sketches a molecule
//top-menu: Chem | Marvin Sketch
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