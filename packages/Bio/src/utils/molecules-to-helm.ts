import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {ToHelmResult} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getMacromoleculeAssembler, setAssemblerLibrary} from '@datagrok-libraries/bio/src/monomer-works/assembler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

/* Molecules to HELM (Bio | Transform | Molecules to HELM..., Molecule to HELM Single, Molecules to HELM Column): the
 * molecules converted by Chem's macromolecule assembler (`IMacromoleculeAssembler.toHelm`: Crux, a column on Chem's
 * workers, a few molecules on the main thread) with the effective monomer library, under the key and JSON To Atomic
 * Level uses (library-json.ts), so a worker compiles the library once for both directions. It replaces the Python
 * script `scripts/mol-to-helm.py`. */

/** The name of the result column (an unused one in the table: `regenerated sequences (2)`, ...), as the script's. */
export const REGENERATED_COL_NAME = 'regenerated sequences';
/** From this many molecules on, the column conversion shows a progress indicator. */
const PROGRESS_FROM_ROWS = 200;

/** The molecules (MOL blocks, SMILES or CXSMILES; null or '' an empty row) converted to HELM with `monomerLib`: per
 * row the HELM ('' for an empty molecule, null where the row has an error), the error and the counts. A molecule that
 * does not convert gets its error; nothing is thrown for it.
 * @param {(string | null)[]} molecules
 * @param {IMonomerLibBase} monomerLib the library to match with
 * @param {(done: number, total: number) => void} onProgress told the molecules done as a column converts
 * @return {Promise<ToHelmResult>} */
export async function convertMoleculesToHelm(
  molecules: (string | null)[], monomerLib: IMonomerLibBase, onProgress?: (done: number, total: number) => void,
): Promise<ToHelmResult> {
  const assembler = await getMacromoleculeAssembler();
  if (typeof assembler.toHelm !== 'function') {
    throw new Error('Molecules to HELM needs a Chem package with the molecule to HELM converter ' +
      '(the toHelm of Chem:getMacromoleculeAssembler): update Chem');
  }
  const libKey = await setAssemblerLibrary(assembler, monomerLib);
  return await assembler.toHelm(molecules, libKey, {onProgress});
}

/** The HELM column of a molecule column, as the script's was: `regenerated sequences` (an unused name), units helm,
 * semantic type Macromolecule, the HELM renderer; '' where a molecule is empty or does not convert. The column is
 * added to `table`. With `notify`, a long column shows its progress and a warning says how many molecules failed.
 * @param {DG.DataFrame} table the molecules' table; gets the column
 * @param {DG.Column<string>} molecules a Molecule column
 * @param {IMonomerLibBase} monomerLib the library to match with
 * @param {boolean} notify the progress indicator and the warning
 * @return {Promise<{helmCol: DG.Column<string>, errors: (string | null)[]}>} the column and the per-row errors */
export async function addHelmColumnOfMolecules(
  table: DG.DataFrame, molecules: DG.Column<string>, monomerLib: IMonomerLibBase, notify: boolean = true,
): Promise<{ helmCol: DG.Column<string>, errors: (string | null)[] }> {
  const values: (string | null)[] = Array.from({length: molecules.length}, (_, i) => molecules.get(i) ?? null);
  const pi = notify && values.length >= PROGRESS_FROM_ROWS ?
    DG.TaskBarProgressIndicator.create('Molecules to HELM...') : null;
  let res: ToHelmResult;
  try {
    res = await convertMoleculesToHelm(values, monomerLib, pi ?
      (done, total) => pi.update(100 * done / Math.max(1, total), `Molecules to HELM: ${done} of ${total}`) :
      undefined);
  } finally {
    pi?.close();
  }
  const name = table.columns.getUnusedName(REGENERATED_COL_NAME);
  const helmCol = DG.Column.fromStrings(name, res.helms.map((h) => h ?? '')) as DG.Column<string>;
  helmCol.meta.units = NOTATION.HELM;
  helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
  helmCol.setTag('cell.renderer', 'helm');
  table.columns.add(helmCol);
  const failed = res.errors.filter((e) => !!e).length;
  if (notify && failed > 0)
    grok.shell.warning(`Molecules to HELM: ${failed} of ${values.length} molecules did not convert`);
  return {helmCol, errors: res.errors};
}
