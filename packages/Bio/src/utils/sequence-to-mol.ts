import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {ISeqHelper, ToAtomicLevelRes} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {buildMonomerHoverLink} from '@datagrok-libraries/bio/src/monomer-works/monomer-hover';
import {assembleSeqColumn} from '@datagrok-libraries/bio/src/monomer-works/assembler';

import {checkInputColumn} from './check-input-column';
import {getMolColName} from '@datagrok-libraries/bio/src/monomer-works/utils';

/** Converts a macromolecule column (any notation) to a molecule column `molfile(<seq>)` added to the table, with the
 * monomer hover link. Both modes run on Chem's macromolecule assembler: `nonlinear` picks the depiction, else the
 * linear layout (Datagrok's poses), which falls back per row to the non-linear one for rows that are not plain
 * chains. The sequence column is read as it is (no notation conversion, no temporary columns).
 * @param {DG.DataFrame} df the table of the sequence column; gets the molecule column
 * @param {DG.Column} macroMolecule the sequence column
 * @param {boolean} nonlinear the non-linear layout for every row
 * @param {boolean} highlight the `.sequence-src-highlight-monomers` tag: highlight the monomers of every molecule
 * @param {IMonomerLibBase} monomerLib the library to resolve with (the column's override library included)
 * @param {ISeqHelper} seqHelper
 * @param {boolean} notify the warning balloon when many rows fail; off for the one-sequence panels
 * @return {Promise<ToAtomicLevelRes>} the molecule column and one warning per row that could not be built */
export async function sequenceToMolfile(
  df: DG.DataFrame, macroMolecule: DG.Column, nonlinear: boolean, highlight: boolean,
  monomerLib: IMonomerLibBase, seqHelper: ISeqHelper, notify: boolean = true,
): Promise<ToAtomicLevelRes> {
  if (!checkInputColumn(macroMolecule, 'To Atomic Level', seqHelper)[0])
    return {molCol: null, warnings: ['Column is not suitable']};

  const res = await assembleSeqColumn(macroMolecule, monomerLib, seqHelper, nonlinear ? 'nonlinear' : 'linear',
    notify);

  if (res.molCol) {
    const molColName = getMolColName(df, macroMolecule.name);
    res.molCol.name = molColName;
    df.columns.add(res.molCol, true);

    await buildMonomerHoverLink(macroMolecule, res.molCol, monomerLib, seqHelper, null, nonlinear);
    res.molCol.setTag(ChemTags.SEQUENCE_SRC_HL_MONOMERS, String(highlight));
    await grok.data.detectSemanticTypes(df);
  }
  return res;
}
