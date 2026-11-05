import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {IMonomerLibBase} from '../types/monomer-library';
import {ISeqHelper, ToAtomicLevelRes} from '../utils/seq-helper';
import {assembleSeqColumn} from './assembler';

/* The linear "To Atomic Level" conversion: {@link _toAtomicLevel} calls Chem's macromolecule assembler (Crux, PLAN
 * D2/D5) with the linear layout. */

/** Convert Macromolecule column into Molecule column storing molfile V3000 with the help of a monomer library:
 * the linear layout (Datagrok's zig-zag / strand poses), and the non-linear one for rows that are not plain chains
 * (cycles, cross-links, CHEM monomers), decided per row by the assembler. Rows that cannot be built give an empty
 * cell and a warning each.
 * @param {DG.DataFrame} df - unused (the result column gets an unused name in seqCol's table)
 * @param {DG.Column} seqCol - Column containing the macromolecule sequence
 * @param {IMonomerLib} monomerLib - Monomer library
 * @param {ISeqHelper} seqHelper - Sequence helper
 * @param {RDModule} _rdKitModule - deprecated and unused: RDKit no longer takes part in the conversion
 * @return {Promise<ToAtomicLevelRes>} the molecule column (not added to the table) and per-row warnings */
export async function _toAtomicLevel(
  df: DG.DataFrame, seqCol: DG.Column<string>,
  monomerLib: IMonomerLibBase, seqHelper: ISeqHelper, _rdKitModule?: RDModule | null,
): Promise<ToAtomicLevelRes> {
  return assembleSeqColumn(seqCol, monomerLib, seqHelper, 'linear');
}
