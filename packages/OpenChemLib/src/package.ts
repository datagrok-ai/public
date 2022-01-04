/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OpenChemLibSketcher} from "./ocl-sketcher";
import {_importSdf} from "./sdf-importer";
import { OCLCellRenderer } from './ocl-cell-renderer';

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

//name: importSdfs
//description: Opens SDF file
//tags: file-handler
//meta.ext: sdf
//input: list bytes
//output: list tables
export function importSdf(bytes: Uint8Array) {
  return _importSdf(Uint8Array.from(bytes));
}

//name: oclCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//meta-cell-renderer-sem-type: Molecule
//output: grid_cell_renderer result
export async function oclCellRenderer() {
  return new OCLCellRenderer();
}
