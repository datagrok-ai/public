import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {previewNglUI, viewNglUI} from './ngl-ui';
import {previewMolstarUI, viewMolstarUI} from './molstar-viewer/utils';
import {BuiltInTrajectoryFormat} from 'molstar/lib/mol-plugin-state/formats/trajectory';

// -- View --

/** View Biostructure (PDB) for file handler
 * @param {string} content - string content of the file
 * @param {BuiltInTrajectoryFormat} format - format of the file
 */
export async function viewBiostructure(content: string, format?: BuiltInTrajectoryFormat): Promise<void> {
  await viewMolstarUI(content, undefined, format);
}

export async function viewNgl(content: string): Promise<void> {
  await viewNglUI(content);
}


// -- Preview --

/** Preview Biostructure (PDB) for file (pre) viewer
 * @param  {DG.FileInfo} file - file to preview, supported formats:
 *  "mmcif" | "cifCore" | "pdb" | "pdbqt" | "gro" | "xyz" | "mol" | "sdf" | "mol2"
 * @return {DG.View} - view of the file
 */
export function previewBiostructure(file: DG.FileInfo): DG.View {
  //return nglViewUI(file); // Deprecated functionality with NGL
  const {view, loadingPromise} = previewMolstarUI(file);
  loadingPromise.catch((err) => {
    const errMsg: string = err instanceof Error ? err.message : err.toString();
    grok.shell.error(`Preview file '' error:\n  ${errMsg}`);
  });
  return view;
}

/** Preview file with NGL, from NglViewer package.
 * @returns {DG.View} */
export function previewNgl(file: DG.FileInfo): DG.View {
  const {view, loadingPromise} = previewNglUI(file);
  loadingPromise.catch((err) => {
    const errMsg: string = err instanceof Error ? err.message : err.toString();
    grok.shell.error(`Preview file '' error:\n  ${errMsg}`);
  });
  return view;
}
