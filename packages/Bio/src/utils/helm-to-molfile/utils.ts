import * as DG from 'datagrok-api/dg';

import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {assembleSeqColumn} from '@datagrok-libraries/bio/src/monomer-works/assembler';
import {getMolColName} from '@datagrok-libraries/bio/src/monomer-works/utils';

import {asHelmColumn} from './converter/converter';

import {_package} from '../../package';

/** Translate HELM column into molfile column (non-linear layout, Chem's macromolecule assembler); the column is not
 * added to the dataframe. chiralityEngine has no effect: stereo is always absolute. The rows are assembled on Chem's
 * workers (as SeqHelper.helmToAtomicLevel does): the main thread stays free, and a row too large for memory fails
 * alone instead of switching Crux off for the session. The molecules are the legacy synchronous converter's. */
export async function getMolColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>, _chiralityEngine: boolean = true, monomerLib: IMonomerLibBase,
): Promise<DG.Column<string>> {
  const res = await assembleSeqColumn(asHelmColumn(helmCol), monomerLib, _package.seqHelper, 'nonlinear');
  const molCol = res.molCol!; // asHelmColumn gives a Macromolecule column, which always gets one
  molCol.name = getMolColName(helmCol.dataFrame, helmCol.name);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  return molCol;
}
