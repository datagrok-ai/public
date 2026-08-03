import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

import * as C from './constants';
import {getMacromoleculeColumns} from './ui-utils';

/** Builds a `MacromoleculeDifference` column from two sequence columns by encoding each row as
 *  `seq1#seq2`. The result copies the relevant metadata (notation, separator, notation provider)
 *  from `seqCol1` — callers are responsible for ensuring the two input columns are compatible. */
export function compareSequencesDo(
  seqCol1: DG.Column<string>, seqCol2: DG.Column<string>, resultName?: string,
): DG.Column<string> {
  const rowCount = Math.min(seqCol1.length, seqCol2.length);
  const values: string[] = new Array(rowCount);
  for (let i = 0; i < rowCount; i++) {
    const v1 = seqCol1.get(i) ?? '';
    const v2 = seqCol2.get(i) ?? '';
    values[i] = `${v1}#${v2}`;
  }

  const name = resultName ?? `${seqCol1.name} vs ${seqCol2.name}`;
  const diffCol = DG.Column.fromStrings(name, values);
  diffCol.semType = C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;

  const sep = seqCol1.getTag(bioTAGS.separator);
  if (sep != null) diffCol.tags[bioTAGS.separator] = sep;
  const units = seqCol1.getTag(DG.TAGS.UNITS);
  if (units != null) diffCol.tags[DG.TAGS.UNITS] = units;
  const aligned = seqCol1.getTag(bioTAGS.aligned);
  if (aligned != null) diffCol.tags[bioTAGS.aligned] = aligned;
  const alphabet = seqCol1.getTag(bioTAGS.alphabet);
  if (alphabet != null) diffCol.tags[bioTAGS.alphabet] = alphabet;

  diffCol.tags[DG.TAGS.CELL_RENDERER] = C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;
  diffCol.temp[SeqTemps.notationProvider] = seqCol1.temp[SeqTemps.notationProvider];

  return diffCol;
}

/** Validates that two macromolecule columns are compatible for pair-wise difference rendering:
 *  identical notation/alignment/alphabet and (for separator notation) matching separator. */
function checkColumnsCompatible(c1: DG.Column, c2: DG.Column): string | null {
  if (c1 === c2) return 'Please choose two distinct columns';
  if (c1.semType !== DG.SEMTYPE.MACROMOLECULE || c2.semType !== DG.SEMTYPE.MACROMOLECULE)
    return 'Both columns must be Macromolecule semantic type';
  if (c1.getTag(DG.TAGS.UNITS) !== c2.getTag(DG.TAGS.UNITS))
    return 'Columns must use the same notation (units)';
  if ((c1.getTag(bioTAGS.separator) ?? '') !== (c2.getTag(bioTAGS.separator) ?? ''))
    return 'Columns must use the same separator';
  if ((c1.getTag(bioTAGS.alphabet) ?? '') !== (c2.getTag(bioTAGS.alphabet) ?? ''))
    return 'Columns must use the same alphabet';
  return null;
}

/** Top-menu entry point for `Bio | Analyze | Compare sequences`. Opens a simple dialog to pick
 *  two Macromolecule columns from the current table and appends a `MacromoleculeDifference`
 *  column to the dataframe. */
export function compareSequencesUI(): void {
  const tv = grok.shell.tv;
  if (!tv || !tv.dataFrame) {
    grok.shell.error('No active table');
    return;
  }
  const df = tv.dataFrame;
  const cols: DG.Column[] = getMacromoleculeColumns();
  if (cols.length < 2) {
    grok.shell.error('Current table needs at least two Macromolecule columns');
    return;
  }

  const names = cols.map((c) => c.name);
  const col1Input = ui.input.choice('Sequence column 1', {value: names[0], items: names});
  const col2Input = ui.input.choice('Sequence column 2', {value: names[1], items: names});
  const resultNameInput = ui.input.string('Result column name', {value: ''});
  resultNameInput.setTooltip('Leave empty to auto-generate from the chosen columns');

  ui.dialog({title: 'Compare Sequences'})
    .add(col1Input)
    .add(col2Input)
    .add(resultNameInput)
    .onOK(() => {
      const c1 = df.col(col1Input.value ?? '');
      const c2 = df.col(col2Input.value ?? '');
      if (!c1 || !c2) {
        grok.shell.error('Could not resolve chosen columns');
        return;
      }
      const err = checkColumnsCompatible(c1, c2);
      if (err != null) {
        grok.shell.error(err);
        return;
      }
      const desired = (resultNameInput.value ?? '').trim() || `${c1.name} vs ${c2.name}`;
      const name = df.columns.getUnusedName(desired);
      const diffCol = compareSequencesDo(
        c1 as DG.Column<string>, c2 as DG.Column<string>, name);
      df.columns.add(diffCol);
    })
    .show();
}
