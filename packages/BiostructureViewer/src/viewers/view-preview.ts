import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {nglViewUI} from './ngl-ui';
import {previewMolstarUI, viewMolstarUI} from './molstar-viewer';

/** View Biostructure (PDB) for file handler */
export async function viewBiostructure(content: string): Promise<void> {
  await viewMolstarUI(content);
}

/** Preview Biostructure (PDB) for file (pre) viewer */
export function previewBiostructure(file: DG.FileInfo): DG.View {
  //return nglViewUI(file); // Deprecated functionality with NGL
  return previewMolstarUI(file);
}
