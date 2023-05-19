import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {nglViewUI} from './ngl-ui';
import {previewMolstarUI, viewMolstarUI} from './molstar-viewer';
import {BuiltInTrajectoryFormat} from 'molstar/lib/mol-plugin-state/formats/trajectory';

/** View Biostructure (PDB) for file handler */
export async function viewBiostructure(content: string, format?: BuiltInTrajectoryFormat): Promise<void> {
  await viewMolstarUI(content, undefined, format);
}

/** Preview Biostructure (PDB) for file (pre) viewer
 * @param file - file to preview, supported formats:
 *  "mmcif" | "cifCore" | "pdb" | "pdbqt" | "gro" | "xyz" | "mol" | "sdf" | "mol2"
*/
export function previewBiostructure(file: DG.FileInfo): DG.View {
  //return nglViewUI(file); // Deprecated functionality with NGL
  return previewMolstarUI(file).view;
}
