/** Profile HMM logos, drawn by Bio's WebLogo (its Profile property): this
 * module supplies the residue heights of every match position exactly as C
 * `hmmlogo` computes them (relative entropy in bits; tools module,
 * hmmer/tools.ts) through a `webLogoProfile` function (package.ts). A profile
 * is one built in the table (Build Profile HMM), a Pfam family of the built-in
 * library or a Pfam accession, or an ANARCI germline model (`germline:human_H`). */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {WebLogoProfile} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {HmmerPool} from './pool';
import {fetchPfam, loadBuiltIn, loadGermlines, loadText, type LoadedModels} from './models';
import {HMMER_NOTE, noteDiv} from './notice';

const LETTERS: Record<string, string> = {amino: 'ACDEFGHIKLMNPQRSTVWY', dna: 'ACGT', rna: 'ACGU'};
const ALPHABETS: Record<string, string> = {amino: 'PT', dna: 'DNA', rna: 'RNA'};
/** Table tags of a profile built in the table: its HMMER text, and the column it was built from. */
const BUILT_TAG = '.hmm.profile.';
const BUILT_COLUMN_TAG = '.hmm.profileColumn.';
/** Tag of an hmmalign column naming its profile, "name accession" (align.ts). */
export const PROFILE_TAG = '.hmmProfile';

/** Keep a profile built from `column` with the table (saved with projects). */
export function storeBuiltProfile(table: DG.DataFrame, name: string, text: string, column: DG.Column): void {
  table.setTag(BUILT_TAG + name, text);
  table.setTag(BUILT_COLUMN_TAG + name, column.name);
}

/** Names of the profiles built in the table. */
export function builtProfiles(table: DG.DataFrame): string[] {
  return Array.from(table.tags.keys()).filter((k) => k.startsWith(BUILT_TAG)).map((k) => k.slice(BUILT_TAG.length));
}

/** Whether `column` is an hmmalign alignment to the profile named `spec` (name or accession). */
function alignedTo(column: DG.Column, spec: string): boolean {
  const [name, accession] = (column.getTag(PROFILE_TAG) ?? '').toLowerCase().split(' ');
  const lower = spec.toLowerCase();
  return !!name && (name === lower || !!accession && accession.split('.')[0] === lower);
}

/** A profile spec as model data and index. */
export async function profileModels(pool: HmmerPool, pkg: DG.Package, profile: string,
  table: DG.DataFrame | null = null): Promise<{models: LoadedModels; index: number}> {
  const spec = profile.trim();
  if (!spec) throw new Error('Choose a profile: a Pfam family (V-set), accession (PF07686) or germline:human_H');
  const built = table?.getTag(BUILT_TAG + spec);
  if (built) return {models: await loadText(pool, built), index: 0};
  // ANARCI germline models are named species_chain (human_H); `germline:` forces them.
  const germline = /^(germline:)?([a-z]+_[HKLABGD])$/.exec(spec);
  if (germline) {
    const germlines = await loadGermlines(pool, pkg);
    const index = germlines.info.models.findIndex((m) => m.name === germline[2]);
    if (index >= 0 || germline[1]) {
      if (index < 0) throw new Error(`No ANARCI germline model ${germline[2]}`);
      return {models: germlines, index};
    }
  }
  const library = await loadBuiltIn(pool, pkg);
  const lower = spec.toLowerCase();
  const index = library.info.models.findIndex((m) =>
    m.name.toLowerCase() === lower || m.accession.split('.')[0].toLowerCase() === lower);
  if (index >= 0) return {models: library, index};
  const accession = /^PF\d{5}$/i.exec(spec)?.[0];
  if (!accession) throw new Error(`${spec} is not in the built-in library; give its Pfam accession (PF00000)`);
  return {models: await fetchPfam(pool, [accession.toUpperCase()]), index: 0};
}

/** Alignment column (1-based; 0 when absent) of every match position: the MAP annotation of HMMER3 text. */
export function mapColumns(text: string): number[] {
  const body = text.slice(Math.max(0, text.search(/^HMM\s/m)));
  return Array.from(body.matchAll(/^\s+\d+(?:\s+\S+){20}\s+(\d+|-)\s/gm), (m) => m[1] === '-' ? 0 : Number(m[1]));
}

/** Position of `column` aligned to each match position of the profile, or undefined
 * when the column is not aligned to it: an hmmalign column of this profile
 * (positions named by match position) or the alignment it was built from (MAP). */
function columnPositions(table: DG.DataFrame | null, column: string, spec: string, length: number):
  number[] | undefined {
  const col = column && table ? table.col(column) : null;
  if (!col || !table) return undefined;
  const built = table.getTag(BUILT_TAG + spec);
  if (built) {
    if (table.getTag(BUILT_COLUMN_TAG + spec) !== column) return undefined;
    const map = mapColumns(built);
    return map.length === length ? map.map((c) => c - 1) : undefined;
  }
  if (!alignedTo(col, spec)) return undefined;
  const names = (col.getTag('.positionNames') ?? '').split(',').map((s) => s.trim());
  return Array.from({length}, (_, i) => names.indexOf(String(i + 1)));
}

/** The WebLogo profile of `profile`: hmmlogo heights (bits, `--height_relent_all`)
 * of every match position, with the model's consensus and insert values. */
export async function webLogoProfile(pool: HmmerPool, pkg: DG.Package, table: DG.DataFrame | null, column: string,
  profile: string): Promise<WebLogoProfile> {
  const {models, index} = await profileModels(pool, pkg, profile, table);
  const model = models.info.models.length === 1 ? models.bytes : await pool.extract(models.bytes, index);
  const [logo] = await pool.logos(model, 'relent_all');
  const info = models.info.models[index];
  const letters = LETTERS[models.info.alphabet] ?? LETTERS.amino;
  return {
    title: `${info.name}${info.accession ? ` (${info.accession})` : ''}, hmmlogo heights (Rusty-HMMER)`,
    alphabet: ALPHABETS[models.info.alphabet] ?? 'PT',
    units: 'bits',
    positions: logo.heights.map((row, i) => ({
      name: String(i + 1),
      heights: Object.fromEntries(Array.from(row, (h, r) => [letters[r], h])),
      info: {
        'Consensus': info.consensus[i + 1] ?? '',
        'Relative entropy': `${logo.relativeEntropy[i].toFixed(3)} bits`,
        'Insert probability': logo.insertProbability[i].toFixed(3),
        'Expected insert length': logo.insertLength[i].toFixed(3),
        'Occupancy': logo.occupancy[i].toFixed(3),
      },
    })),
    columnPositions: columnPositions(table, column, profile.trim(), logo.heights.length),
  };
}

/** The column of `table` aligned to the profile `spec`, if any (see {@link columnPositions}). */
export function alignedColumn(table: DG.DataFrame, spec: string): DG.Column | undefined {
  const builtFrom = table.getTag(BUILT_COLUMN_TAG + spec);
  return builtFrom ? table.col(builtFrom) ?? undefined : table.columns.toList().find((c) => alignedTo(c, spec));
}

/** Add Bio's WebLogo drawing `profile` to `tv`, docked below the grid, for the
 * sequence column `column` (default: the one aligned to the profile, if any);
 * when that column is aligned to the profile, clicking a letter selects the
 * rows that have it there. */
export function addProfileLogo(tv: DG.TableView, profile: string, column?: DG.Column | null): DG.Viewer {
  column ??= alignedColumn(tv.dataFrame, profile);
  const viewer = tv.addViewer('WebLogo', {profile, ...(column ? {sequenceColumnName: column.name} : {})});
  tv.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, `${profile} logo`, 0.3);
  return viewer;
}

const SOURCES = ['Pfam family (built-in)', 'Pfam accession (online)', 'ANARCI germline', 'Built in this table'];
const GERMLINE = /^[a-z]+_[HKLABGD]$/;

/** The profile a column is aligned to, as a dialog source and value: an hmmalign
 * column (Pfam family or germline) or the alignment a profile was built from. */
function columnProfile(table: DG.DataFrame, column: DG.Column | null): [string, string] | null {
  if (!column) return null;
  const built = builtProfiles(table).filter((name) => table.getTag(BUILT_COLUMN_TAG + name) === column.name);
  if (built.length > 0) return [SOURCES[3], built[built.length - 1]];
  const name = column.getTag(PROFILE_TAG)?.split(' ')[0];
  if (!name) return null;
  return GERMLINE.test(name) ? [SOURCES[2], name] : [SOURCES[0], name];
}

/** `Bio | Analyze | Profile HMM Logo...`: choose a sequence column and a profile, and add its WebLogo. */
export async function showLogoDialog(pool: HmmerPool, pkg: DG.Package): Promise<void> {
  const tv = grok.shell.tv;
  const df = tv?.dataFrame;
  const seqCols = df?.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE) ?? [];
  if (!tv || !df || seqCols.length === 0) {
    grok.shell.warning('Open a table with a Macromolecule column');
    return;
  }
  const built = builtProfiles(df);
  const sources = SOURCES.filter((s) => s !== SOURCES[3] || built.length > 0);
  // Default column: one aligned to a profile (hmmalign, or a profile built from it), else the first.
  const columnInput = ui.input.column('Sequence', {table: df,
    value: seqCols.find((c) => columnProfile(df, c)) ?? seqCols[0],
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MACROMOLECULE,
    tooltipText: 'The WebLogo\'s sequence column. When it is aligned to the profile (an hmmalign column of it, ' +
      'or the alignment it was built from), clicking a letter selects the rows that have it there.'});
  const sourceInput = ui.input.choice('Profile from', {value: sources[0], items: sources});
  const familyInput = ui.input.choice<string>('Family', {value: '', items: [''], nullable: false,
    tooltipText: 'A Pfam family of the built-in library'});
  const accessionInput = ui.input.string('Accession', {value: 'PF07686',
    tooltipText: 'Any Pfam accession, downloaded from InterPro (e.g. PF07686)'});
  const germlineInput = ui.input.choice<string>('Germline', {value: '', items: [''], nullable: false,
    tooltipText: 'An ANARCI germline model: species and chain (H heavy, K kappa, L lambda, A/B/G/D TCR)'});
  const builtInput = ui.input.choice<string>('Profile', {value: built[built.length - 1] ?? '', items: built,
    tooltipText: 'A profile built from rows of this table (Build Profile HMM)'});
  const valueInputs: Record<string, DG.InputBase> = {[SOURCES[0]]: familyInput, [SOURCES[1]]: accessionInput,
    [SOURCES[2]]: germlineInput, [SOURCES[3]]: builtInput};
  let library: LoadedModels | null = null;
  let germlines: LoadedModels | null = null;
  const names = (models: LoadedModels) => models.info.models.map((m) => m.name);
  const refresh = async () => {
    for (const [source, input] of Object.entries(valueInputs))
      input.root.style.display = source === sourceInput.value ? '' : 'none';
    if (sourceInput.value === SOURCES[0] && !library) {
      library = await loadBuiltIn(pool, pkg);
      familyInput.items = names(library).sort((a, b) => a.localeCompare(b));
      familyInput.value ||= 'V-set';
    }
    if (sourceInput.value === SOURCES[2] && !germlines) {
      germlines = await loadGermlines(pool, pkg);
      germlineInput.items = names(germlines);
      germlineInput.value ||= 'human_H';
    }
  };
  // Picking a column aligned to a profile picks that profile.
  const followColumn = async () => {
    const aligned = columnProfile(df, columnInput.value);
    if (!aligned) return;
    const [source, value] = aligned;
    if (source === SOURCES[0]) familyInput.value = value;
    if (source === SOURCES[2]) germlineInput.value = value;
    if (source === SOURCES[3]) builtInput.value = value;
    sourceInput.value = source;
    await refresh();
    // A family or germline appears in its list once loaded.
    if (source === SOURCES[0]) familyInput.value = value;
    if (source === SOURCES[2]) germlineInput.value = value;
  };
  const report = (e: unknown) => grok.shell.error(`${e instanceof Error ? e.message : e}`);
  sourceInput.onChanged.subscribe(() => refresh().catch(report));
  columnInput.onChanged.subscribe(() => followColumn().catch(() => {}));
  await refresh();
  await followColumn();

  ui.dialog({title: 'Profile HMM Logo'})
    .add(ui.inputs([columnInput, sourceInput, familyInput, accessionInput, germlineInput, builtInput]))
    .add(ui.divText('Drawn by WebLogo with hmmlogo heights; change the profile later in its Profile property.',
      {style: {fontSize: '12px', opacity: '0.7', maxWidth: '360px'}}))
    .add(noteDiv(HMMER_NOTE))
    .onOK(() => {
      const value = valueInputs[sourceInput.value ?? SOURCES[0]]?.value;
      const spec = String(value ?? '').trim();
      if (!spec) {
        grok.shell.warning('Choose a profile');
        return;
      }
      addProfileLogo(tv, spec, columnInput.value);
    })
    .show();
}
