import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { BioStructureViewer } from './biostrucure-viewer';

export const _package = new DG.Package();

//tags: fileViewer, fileViewer-mol, fileViewer-mol2, fileViewer-cif, fileViewer-mcif, fileViewer-mmcif, fileViewer-gro, fileViewer-pdb, fileViewer-ent, fileViewer-pqr, fileViewer-mmtf, fileViewer-mtl, fileViewer-sdf, fileViewer-sd
//input: file file
//output: view v
export function nglStructureViewer(file: any) {
  return nglViewer(file);
}

//tags: fileViewer, fileViewer-ply, fileViewer-obj
//input: file file
//output: view v
export function nglSurfaceViewer(file: any) {
  return nglViewer(file);
}

//tags: fileViewer, fileViewer-prmtop, fileViewer-parm7, fileViewer-psf, fileViewer-top
//input: file file
//output: view v
export function nglTopologyViewer(file: any) {
  return nglViewer(file);
}

//tags: fileViewer, fileViewer-dsn6, fileViewer-brix, fileViewer-cube, fileViewer-cub, fileViewer-dx, fileViewer-dxbin, fileViewer-xplor, fileViewer-cns, fileViewer-mrc, fileViewer-map, fileViewer-ccp4
//input: file file
//output: view v
export function nglDensityViewer(file: any) {
  return nglViewer(file);
}

//name: PDB Viewer
//tags: panel
//input: string pdbId {semType: pdb_id}
//output: widget w
export function nglPdbPanelWidget(pdbId: string) {
  var host = ui.div([], 'd4-ngl-viewer');
  //@ts-ignore
  var stage = new NGL.Stage(host);
  handleResize(host, stage);
  stage.loadFile(`rcsb://${pdbId}`, { defaultRepresentation: true });
  return DG.Widget.fromRoot(host);
}

function handleResize(host: HTMLDivElement, stage: any) {
  let canvas = host.querySelector('canvas');
  function resize() {
    canvas!.width = Math.floor(canvas!.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(canvas!.clientHeight * window.devicePixelRatio);
    stage.handleResize();
  }

  ui.onSizeChanged(host).subscribe((_) => resize());
  resize();
}


/** @returns {DG.View} */
export default function nglViewer(file: any): DG.View {
  let view = DG.View.create();
  var host = ui.div([], 'd4-ngl-viewer');
  //@ts-ignore
  var stage = new NGL.Stage(host);
  handleResize(host, stage);

  function loadBytes(bytes: any) {
    var blob = new Blob([bytes], {type: 'application/octet-binary'});
    stage.loadFile(blob, { defaultRepresentation: true, ext: file.extension} );
  }

  file
    .readAsBytes()
    .then(loadBytes);

  view.append(host);
  return view;
}

//name: Docking
//tags: app
export async function biostructureApp() {

  let pi = DG.TaskBarProgressIndicator.create('Opening BioStructure Viewer');
  let app = new BioStructureViewer();
  await app.init();
  pi.close();
}