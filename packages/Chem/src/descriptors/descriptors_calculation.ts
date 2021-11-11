import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

export function getDescriptors(smiles: DG.Column) {
    return(new DG.Column(smiles));
}