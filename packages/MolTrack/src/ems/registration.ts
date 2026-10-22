import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {domains} from '@datagrok-libraries/ems/src/dg';
import type {DomainTable} from '@datagrok-libraries/ems/src/dg';
import {badge} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, HASH_COLUMNS, IDENTITY_COLUMNS} from './chemistry';
import type {CompoundIdentity} from './chemistry';
import {MolTrackTables} from './tables';

const STRUCTURE_NAMES = ['smiles', 'structure', 'molecule', 'molfile', 'molblock', 'original_structure'];

const LOOKUP_CHUNK = 500;

const PREVIEW_ROWS = 50;

export interface PreparedRegistration {
  ready: DG.DataFrame;
  /** '' where the structure cannot be read */
  hashes: string[];
  structure: DG.Column;
}

type Verdict = 'New' | 'Existing' | 'Repeated' | 'Error';

export class BulkRegistration {
  static readonly VERDICT = 'Verdict';
  static readonly RESULT = 'Registration';
  static readonly GROK_ID = 'Grok ID';
  static readonly NEW = 'New compound';
  private static readonly COLORS: Record<Verdict, number> = {New: DG.Color.success, Existing: DG.Color.blue,
    Repeated: DG.Color.gray, Error: DG.Color.failure};
  private static readonly BADGES: Record<Verdict, BadgeVariant> = {New: 'accent', Existing: 'default',
    Repeated: 'default', Error: 'error'};

  // the import maps by caption: `structure` (Parent Structure) and `original_structure` (Structure) stay apart
  static label(compounds: DomainTable<CompoundRow>, column: string): string {
    return compounds.properties.find((p) => p.name === column)?.friendlyName ?? column;
  }

  static async prepare(df: DG.DataFrame, compounds: DomainTable<CompoundRow>): Promise<PreparedRegistration> {
    const column = await BulkRegistration._structureColumn(df);
    const rule = MolTrackTables.matchRule();
    const identities: (CompoundIdentity | null)[] = [];
    for (let i = 0; i < df.rowCount; i++) {
      const input = column.isNone(i) ? '' : String(column.get(i));
      identities.push(input === '' ? null : await Chemistry.identify(input, rule).catch(() => null));
    }
    const hashes = identities.map((identity) => identity?.hash_mol ?? '');
    const registered = await BulkRegistration._registered([...new Set(hashes.filter((h) => h !== ''))], rule);
    const first = new Map<string, number>();
    const marks = hashes.map((hash, i): [Verdict, string] => {
      if (hash === '')
        return ['Error', column.isNone(i) ? 'No structure' : 'Cannot read the structure'];
      const molregno = registered.get(hash);
      if (molregno !== undefined)
        return ['Existing', `Already registered as ${MolTrackTables.corporateId(compounds, 'molregno', molregno)}`];
      if (first.has(hash))
        return ['Repeated', `Same compound as row ${first.get(hash)! + 1}`];
      first.set(hash, i);
      return ['New', BulkRegistration.NEW];
    });
    for (const name of [BulkRegistration.VERDICT, BulkRegistration.RESULT]) {
      if (df.col(name) !== null)
        df.columns.remove(name);
    }
    df.columns.addNewString(BulkRegistration.VERDICT).init((i) => marks[i][0])
      .meta.colors.setCategorical(BulkRegistration.COLORS);
    df.columns.addNewString(BulkRegistration.RESULT).init((i) => marks[i][1]);
    const results = marks.map(([, message]) => message);
    const keep = DG.BitSet.create(df.rowCount, (i) => results[i] === BulkRegistration.NEW);
    const ready = df.clone(keep);
    const kept = [...first.values()];
    const label = (name: string) => BulkRegistration.label(compounds, name);
    for (const name of [column.name, BulkRegistration.VERDICT, BulkRegistration.RESULT, BulkRegistration.GROK_ID,
      'original_structure', ...IDENTITY_COLUMNS].map(label)) {
      if (ready.col(name) !== null)
        ready.columns.remove(name);
    }
    ready.columns.addNewString(label('original_structure')).init((i) => String(column.get(kept[i])));
    for (const name of IDENTITY_COLUMNS)
      ready.columns.addNewString(label(name)).init((i) => String(identities[kept[i]]![name]));
    ready.name = `${df.name} - new compounds`;
    return {ready, hashes, structure: column};
  }

  static run(compounds: DomainTable<CompoundRow>): Promise<boolean> {
    const data = ui.input.table('Data', {value: grok.shell.t ?? grok.shell.tables[0]});
    const summary = ui.divText('');
    const preview = ui.div([], 'moltrack-ems-bulk-preview');
    let prepared: Promise<PreparedRegistration> | null = null;
    const check = () => {
      const df = data.value;
      preview.replaceChildren();
      summary.textContent = df ? 'Checking the structures…' : '';
      prepared = df ? BulkRegistration.prepare(df, compounds) : null;
      prepared?.then((p) => {
        if (data.value !== df)
          return;
        summary.textContent = BulkRegistration.summary(df!);
        preview.append(BulkRegistration.preview(df!, p.structure));
      }, (e) => summary.textContent = e instanceof Error ? e.message : String(e));
    };
    data.onChanged.subscribe(check);
    check();
    return new Promise((resolve) => ui.dialog('Register compounds').add(data).add(summary).add(preview)
      .onOK(() => {
        const df = data.value;
        if (df === null || df === undefined || prepared === null)
          return resolve(false);
        prepared.then((p) => BulkRegistration._register(df, p, compounds)).then(resolve, (e) => {
          grok.shell.error(`Bulk registration failed: ${e instanceof Error ? e.message : e}`);
          resolve(false);
        });
      })
      .onCancel(() => resolve(false))
      .show({resizable: true, width: 760, height: 560}));
  }

  static preview(df: DG.DataFrame, structure: DG.Column): HTMLElement {
    const verdicts = df.col(BulkRegistration.VERDICT)!;
    const results = df.col(BulkRegistration.RESULT)!;
    const name = df.columns.toList().find((c) => c !== structure && /name/i.test(c.name));
    const rows = Array.from({length: Math.min(df.rowCount, PREVIEW_ROWS)}, (_, i) => i);
    return ui.table(rows, (i) => {
      const verdict = verdicts.get(i) as Verdict;
      const drawn = verdict === 'Error' ? ui.divText(structure.isNone(i) ? '' : String(structure.get(i))) :
        grok.chem.drawMolecule(String(structure.get(i)), 120, 50);
      return [String(i + 1), drawn, name?.isNone(i) === false ? String(name.get(i)) : '',
        badge(verdict, {variant: BulkRegistration.BADGES[verdict]}), results.get(i)];
    }, ['#', 'Structure', name?.name ?? '', 'Becomes', '']);
  }

  static summary(df: DG.DataFrame): string {
    const verdicts = df.col(BulkRegistration.VERDICT)!.toList() as Verdict[];
    const count = (v: Verdict) => verdicts.filter((x) => x === v).length;
    const parts = ([[count('New'), 'new'], [count('Existing'), 'already registered'],
      [count('Repeated'), 'repeated in the table'], [count('Error'), 'cannot be read']] as [number, string][])
      .filter(([n]) => n > 0).map(([n, what]) => `${n} ${what}`);
    return `${df.rowCount} rows: ${parts.join(', ')}. Only the new ones are registered.`;
  }

  private static async _register(df: DG.DataFrame, prepared: PreparedRegistration,
    compounds: DomainTable<CompoundRow>): Promise<boolean> {
    if (prepared.ready.rowCount === 0 || await domains.import(compounds, {source: prepared.ready}) === null)
      return false;
    const registered = await BulkRegistration._registered([...new Set(prepared.hashes.filter((h) => h !== ''))],
      MolTrackTables.matchRule());
    if (df.col(BulkRegistration.GROK_ID) !== null)
      df.columns.remove(BulkRegistration.GROK_ID);
    df.columns.addNewString(BulkRegistration.GROK_ID).init((i) => MolTrackTables.corporateId(compounds, 'molregno',
      registered.get(prepared.hashes[i])));
    grok.shell.info(`Grok IDs are in the ${BulkRegistration.GROK_ID} column of ${df.name}`);
    return true;
  }

  private static async _structureColumn(df: DG.DataFrame): Promise<DG.Column> {
    await grok.data.detectSemanticTypes(df);
    const column = df.columns.bySemType(DG.SEMTYPE.MOLECULE) ??
      df.columns.toList().find((c) => STRUCTURE_NAMES.includes(c.name.toLowerCase()));
    if (column === undefined || column === null)
      throw new Error('The file has no structure column');
    return column;
  }

  private static async _registered(hashes: string[], rule: keyof typeof HASH_COLUMNS): Promise<Map<string, number>> {
    const column = HASH_COLUMNS[rule];
    const found = new Map<string, number>();
    for (let at = 0; at < hashes.length; at += LOOKUP_CHUNK) {
      const chunk = hashes.slice(at, at + LOOKUP_CHUNK);
      const rows = await moltrackEmsDb.compounds.query({filter: DG.cond(column, '=', chunk),
        columns: [column, 'molregno'], sort: 'molregno'});
      for (const row of rows) {
        if (!found.has(row[column]!))
          found.set(row[column]!, row.molregno);
      }
    }
    return found;
  }
}
