/** Build a profile HMM from aligned rows of a table (`hmmbuild`, exact C
 * model text from the tools module), kept with the table; then search a
 * sequence column with it (`hmmsearch` of that text, as C would read the file)
 * and annotate the hits; optionally show its logo (Bio's WebLogo) and
 * download the `.hmm` file. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {BuiltProfile} from './hmmer/tools.ts';
import {HmmerPool} from './pool';
import {alignedRows} from './sequences';
import {loadText, type LoadedModels} from './models';
import {applySearch, searchColumn} from './search';
import {HMMER_NOTE, withNotes} from './notice';
import {taskProgress} from './progress';

/** hmmbuild model names are single words. */
export function modelName(text: string): string {
  return text.trim().replace(/\s+/g, '_') || 'profile';
}

function alphabetOf(column: DG.Column): 'amino' | 'dna' | 'rna' {
  const alphabet = column.getTag('alphabet');
  return alphabet === 'DNA' ? 'dna' : alphabet === 'RNA' ? 'rna' : 'amino';
}

/** hmmbuild of `rows` of the aligned `column`. */
export async function buildProfile(pool: HmmerPool, column: DG.Column<string>, rows: number[], name: string):
  Promise<BuiltProfile & {models: LoadedModels}> {
  if (rows.length === 0) throw new Error('No rows to build the profile from');
  const aligned = await alignedRows(column, rows);
  const lengths = new Set(aligned.map((r) => r.length));
  if (lengths.size > 1 || aligned[0].length === 0)
    throw new Error(`${column.name} is not aligned (rows of different lengths): align it first (Bio | Analyze | MSA)`);
  const built = await pool.build(aligned.map((a, i) => ({name: `row${rows[i] + 1}`, aligned: a})), modelName(name),
    alphabetOf(column));
  return {...built, models: await loadText(pool, built.text)};
}

/** `Bio | Search | Build Profile HMM...` */
export async function showBuildProfileDialog(pool: HmmerPool): Promise<void> {
  const df = grok.shell.tv?.dataFrame;
  const seqCols = df?.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE) ?? [];
  if (!df || seqCols.length === 0) {
    grok.shell.warning('Open a table with an aligned Macromolecule column');
    return;
  }
  // A column is aligned when every non-empty row has the same length (FASTA; MSA tags are not
  // reliable, e.g. on unaligned sample data), else the MSA tag (other notations).
  const equalLengths = (c: DG.Column) => c.meta.units === 'fasta' &&
    new Set((c.toList() as (string | null)[]).filter((v) => v).map((v) => v!.length)).size === 1;
  const isAligned = (c: DG.Column) => equalLengths(c) ||
    c.meta.units !== 'fasta' && (c.getTag('aligned') ?? '').includes('MSA');
  const alignedInput = ui.input.column('Alignment', {table: df, value: seqCols.find(isAligned) ?? seqCols[0],
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MACROMOLECULE,
    tooltipText: withNotes('An aligned column (e.g. from Bio | Analyze | MSA); its rows build the profile ' +
      '(hmmbuild).', HMMER_NOTE)});
  const rowChoices = ['Selected rows', 'Filtered rows', 'All rows'];
  const rowsInput = ui.input.choice('Rows', {value: df.selection.trueCount > 0 ? rowChoices[0] : rowChoices[2],
    items: rowChoices});
  const nameInput = ui.input.string('Name', {value: `${alignedInput.value?.name ?? 'column'}_profile`});
  const searchInput = ui.input.bool('Search the table', {value: true,
    tooltipText: 'Search a sequence column with the new profile (hmmsearch) and annotate the hits'});
  const targetInput = ui.input.column('Search in', {table: df, value: seqCols.find((c) => !isAligned(c)) ?? seqCols[0],
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MACROMOLECULE});
  const evalueInput = ui.input.float('E-value ≤', {value: 10});
  const logoInput = ui.input.bool('Show logo', {value: true,
    tooltipText: 'Add a WebLogo of the profile (hmmlogo heights); clicking a letter selects the rows of the ' +
      'alignment that have it at that position'});
  const downloadInput = ui.input.bool('Download .hmm', {value: false});
  searchInput.onChanged.subscribe(() => {
    targetInput.root.style.display = evalueInput.root.style.display = searchInput.value ? '' : 'none';
  });

  ui.dialog({title: 'Build Profile HMM'})
    .add(ui.inputs([alignedInput, rowsInput, nameInput, searchInput, targetInput, evalueInput, logoInput,
      downloadInput]))
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('Building profile HMM...');
      try {
        const column = alignedInput.value as DG.Column<string>;
        const mask = rowsInput.value === rowChoices[0] ? df.selection :
          rowsInput.value === rowChoices[1] ? df.filter : null;
        const rows = mask ? Array.from(mask.getSelectedIndexes()) :
          Array.from({length: df.rowCount}, (_, i) => i);
        const built = await buildProfile(pool, column, rows, nameInput.value ?? 'profile');
        const model = built.models.info.models[0];
        const {addProfileLogo, storeBuiltProfile} = await import('./logo');
        storeBuiltProfile(df, model.name, built.text, column);
        if (downloadInput.value) DG.Utils.download(`${model.name}.hmm`, built.text);
        if (searchInput.value) {
          const target = targetInput.value as DG.Column<string>;
          const result = await searchColumn(pool, built.models, 0, target, {evalue: evalueInput.value ?? 10}, false,
            undefined, taskProgress(pi, 'Searching'));
          applySearch(df, target, built.models, 0, result,
            `${model.name} (hmmbuild of ${rows.length} rows of ${column.name})`);
          grok.shell.info(`${model.name}: ${model.length} positions from ${rows.length} rows; ` +
            `${result.reported} sequences found, ${result.included} significant`);
        } else
          grok.shell.info(`${model.name}: ${model.length} positions from ${rows.length} rows`);
        if (logoInput.value) addProfileLogo(grok.shell.tv, model.name);
      } catch (e) {
        grok.shell.error(`Build profile failed: ${e instanceof Error ? e.message : e}`);
      } finally {
        pi.close();
      }
    })
    .show();
}
