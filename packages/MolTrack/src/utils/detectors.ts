import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const GROK_ID_SEM_TYPE = 'Grok ID';
const GROK_BATCH_ID_SEM_TYPE = 'Grok Batch ID';

export function registerSemanticTypes(): void {
  DG.SemanticValue.registerRegExpDetector(GROK_ID_SEM_TYPE, '^DG-\\d{6}$');
  DG.SemanticValue.registerRegExpDetector(GROK_BATCH_ID_SEM_TYPE, '^DGB-\\d{6}$');
}
