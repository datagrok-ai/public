/*
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {BioStructureViewer} from './biostructure-viewer';
import {previewNgl} from './viewers/view-preview';

export const _package = new DG.Package();

//eslint-disable-next-line max-len
//meta.role: fileViewer
//meta.fileViewer: mol,mol2,cif,mcif,mmcif,grok,pdb,ent,pqr,mmtf,mtl,sdf,sd
//input: file file
//output: view v
export function nglStructureViewer(file: any) {
  return previewNgl(file);
}

//meta.role: fileViewer
//meta.fileViewer: ply,obj
//input: file file
//output: view v
export function nglSurfaceViewer(file: any) {
  return previewNgl(file);
}

//meta.role: fileViewer
//meta.fileViewer: prmtop,parm7,psf,top
//input: file file
//output: view v
export function nglTopologyViewer(file: any) {
  return previewNgl(file);
}

//eslint-disable-next-line max-len
//meta.role: fileViewer
//meta.fileViewer: dsn6,brix,cube,cub,dx,dxbin,xplor,cns,mrc,map,ccp4
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
