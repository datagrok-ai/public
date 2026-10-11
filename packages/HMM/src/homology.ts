/** Homology search over a sequence column ("phmmer against the table"): one
 * query sequence scored against every row with bit scores and E-values, as C
 * phmmer computes them. The calibrated single-sequence query model is built
 * by the tools module (hmmer/tools.ts, `p7_SingleBuilder`) and the column is
 * searched by the engine split across workers (identical to one search). */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {SearchOptions} from './hmmer/engine.ts';
import type {QueryModelOptions} from './hmmer/tools.ts';
import {HmmerPool, type Progress} from './pool';
import {taskProgress} from './progress';
import {cellResidues} from './sequences';
import {loadModels, type LoadedModels} from './models';
import {applySearch, type ColumnSearch, searchColumn} from './search';
import {HMMER_NOTE, withNotes} from './notice';

export const MATRICES = ['BLOSUM62', 'BLOSUM45', 'BLOSUM50', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM120',
  'PAM240'];

/** Query residues as phmmer reads them: letters only, upper case (FASTA text). */
export function queryResidues(text: string): string {
  return text.replace(/^>.*$/m, '').replace(/[\s\-.*]/g, '').toUpperCase();
}

/** phmmer of `query` against `column`; columns are labelled `label`. */
export async function homologySearch(pool: HmmerPool, column: DG.Column<string>, query: string,
  options: SearchOptions & QueryModelOptions = {}, label = 'query', progress?: Progress):
  Promise<{models: LoadedModels; result: ColumnSearch}> {
  const residues = queryResidues(query);
  if (!residues) throw new Error('The query sequence is empty');
  const {popen, pextend, matrix, seed, ...search} = options;
  const bytes = await pool.queryModels([{name: label, residues}], {popen, pextend, matrix, seed});
  const models = await loadModels(pool, `phmmer ${label}`, bytes);
  const result = await searchColumn(pool, models, 0, column, {seed, ...search}, false, label, progress);
  return {models, result};
}

/** `Bio | Search | Homology Search (phmmer)...` */
export async function showHomologySearchDialog(pool: HmmerPool): Promise<void> {
  const df = grok.shell.tv?.dataFrame;
  const seqCols = df?.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE) ?? [];
  if (!df || seqCols.length === 0) {
    grok.shell.warning('Open a table with a Macromolecule column');
    return;
  }
  const columnInput = ui.input.column('Sequence', {table: df, value: seqCols[0],
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MACROMOLECULE});
  const sources = ['Current row', 'Sequence'];
  const sourceInput = ui.input.choice('Query', {value: sources[0], items: sources,
    tooltipText: withNotes('The query protein sequence: a row of the column, or typed (phmmer).', HMMER_NOTE)});
  const rowInput = ui.input.int('Row', {value: Math.max(0, df.currentRowIdx) + 1,
    tooltipText: 'The row whose sequence is the query (1-based)'});
  const textInput = ui.input.textArea('Query sequence', {value: '',
    tooltipText: 'Protein sequence (plain or FASTA)'});
  const evalueInput = ui.input.float('E-value ≤', {value: 10, tooltipText: 'Report sequences up to this E-value (-E)'});
  const matrixInput = ui.input.choice('Matrix', {value: MATRICES[0], items: MATRICES,
    tooltipText: 'Substitution matrix of the query model (--mx)'});
  const popenInput = ui.input.float('Gap open', {value: 0.02, tooltipText: 'Gap open probability (--popen)'});
  const pextendInput = ui.input.float('Gap extend', {value: 0.4, tooltipText: 'Gap extend probability (--pextend)'});
  const refresh = () => {
    rowInput.root.style.display = sourceInput.value === sources[0] ? '' : 'none';
    textInput.root.style.display = sourceInput.value === sources[1] ? '' : 'none';
  };
  sourceInput.onChanged.subscribe(refresh);
  refresh();

  ui.dialog({title: 'Homology Search (phmmer)'})
    .add(ui.inputs([columnInput, sourceInput, rowInput, textInput, evalueInput, matrixInput, popenInput,
      pextendInput]))
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('Homology search...');
      try {
        const column = columnInput.value as DG.Column<string>;
        let query: string;
        let label: string;
        if (sourceInput.value === sources[0]) {
          const row = (rowInput.value ?? 1) - 1;
          if (row < 0 || row >= df.rowCount) throw new Error(`Row ${row + 1} is not in the table`);
          query = (await cellResidues(column, row)).residues;
          label = `phmmer row ${row + 1}`;
        } else {
          query = textInput.value ?? '';
          label = 'phmmer query';
        }
        const {models, result} = await homologySearch(pool, column, query, {evalue: evalueInput.value ?? 10,
          matrix: matrixInput.value ?? 'BLOSUM62', popen: popenInput.value ?? 0.02,
          pextend: pextendInput.value ?? 0.4}, label, taskProgress(pi, 'Searching'));
        applySearch(df, column, models, 0, result, `${label} (phmmer, ${matrixInput.value})`);
        grok.shell.info(`${label}: ${result.reported} sequences found, ${result.included} significant`);
      } catch (e) {
        grok.shell.error(`Homology search failed: ${e instanceof Error ? e.message : e}`);
      } finally {
        pi.close();
      }
    })
    .show();
}
