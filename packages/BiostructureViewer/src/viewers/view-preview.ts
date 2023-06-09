import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {previewMolstarUI, viewMolstarUI} from './molstar-viewer';
import {BuiltInTrajectoryFormat} from 'molstar/lib/mol-plugin-state/formats/trajectory';

/** View Biostructure (PDB) for file handler
 * @param {string} content - string content of the file
 * @param {BuiltInTrajectoryFormat} format - format of the file
 */
export async function viewBiostructure(content: string, format?: BuiltInTrajectoryFormat): Promise<void> {
  await viewMolstarUI(content, undefined, format);
}

/** Preview Biostructure (PDB) for file (pre) viewer
 * @param {DG.FileInfo} file - file to preview
 * @return {DG.View} - view of the file
*/
export function previewBiostructure(file: DG.FileInfo): DG.View {
  //return nglViewUI(file); // Deprecated functionality with NGL
  return previewMolstarUI(file).view;
}
