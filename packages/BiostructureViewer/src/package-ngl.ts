/*
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {BioStructureViewer} from './biostructure-viewer';
import {previewNgl} from './viewers/view-preview';

export const _package = new DG.Package();

//eslint-disable-next-line max-len
//tags: fileViewer, fileViewer-mol, fileViewer-mol2, fileViewer-cif, fileViewer-mcif, fileViewer-mmcif, fileViewer-gro, fileViewer-pdb, fileViewer-ent, fileViewer-pqr, fileViewer-mmtf, fileViewer-mtl, fileViewer-sdf, fileViewer-sd
//input: file file
//output: view v
export function nglStructureViewer(file: any) {
  return previewNgl(file);
}

//tags: fileViewer, fileViewer-ply, fileViewer-obj
//input: file file
//output: view v
export function nglSurfaceViewer(file: any) {
  return previewNgl(file);
}

//tags: fileViewer, fileViewer-prmtop, fileViewer-parm7, fileViewer-psf, fileViewer-top
//input: file file
//output: view v
export function nglTopologyViewer(file: any) {
  return previewNgl(file);
}

//eslint-disable-next-line max-len
//tags: fileViewer, fileViewer-dsn6, fileViewer-brix, fileViewer-cube, fileViewer-cub, fileViewer-dx, fileViewer-dxbin, fileViewer-xplor, fileViewer-cns, fileViewer-mrc, fileViewer-map, fileViewer-ccp4
//input: file file
//output: view v
export function nglDensityViewer(file: any) {
  return previewNgl(file);
}


//name: Docking
//tags: app
export async function biostructureApp() {
  const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure Viewer');
  const app = new BioStructureViewer();
  await app.init();
  pi.close();
}
/**/
