import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {DomainApp, DomainFilters} from '@datagrok-libraries/ems/src/dg';
import type {DomainAppOptions, DomainTable} from '@datagrok-libraries/ems/src/dg';
import type {DomainSource} from '@datagrok-libraries/ems';
import {Filters, Form, Rows, Section, isEmpty} from '@datagrok-libraries/u2';
import {ObjectForm} from '@datagrok-libraries/u2/src/dg/index.js';
import type {Action, RowView} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {Chemistry, IDENTITY_COLUMNS} from './chemistry';
import type {CompoundIdentity} from './chemistry';
import {StructureSearch} from './structure-search';
import {MolTrackTables} from './tables';

interface Verdict {
  input: string;
  pending: boolean;
  problem: string | null;
}

// EMS has no server write hook: a draft's identity is derived in the browser as its structure is drawn
export class CompoundsApp extends DomainApp {
  static readonly FIELDS = ['molregno', 'original_structure', 'common_name', 'cas', 'epa_compound_id', 'usan',
    'iupac_name', 'source', 'source_compound_code', 'mol_logp'];

  private readonly _verdicts = new Map<string, Verdict>();
  private readonly _identity = this.runInScope(() => new Section({title: 'Identity', collapsible: false}));
  private _identityForm: Form | null = null;

  constructor(options: DomainAppOptions) {
    super(options);
    this._identity.root.classList.add('moltrack-ems-identity');
    this.own(DomainFilters.registerOperators(this.table.address, StructureSearch.operators()));
    this.own((this.table as DomainTable<CompoundRow>).validators.add('original_structure', (_, row) =>
      this._verdicts.get(row.id)?.problem ?? null));
    this.listSource.hosted = (condition, signal) => StructureSearch.matching(condition, signal);
    // a saved compound shows the identity stored with it; only a draft's is derived
    this.effect(() => {
      const source = this.entitySource.value;
      if (source === null || source.isDraft)
        return;
      const row = source.currentRow.value as RowView<CompoundRow> | null;
      if (row === null)
        return;
      this._showIdentity(row, null, row.original_structure);
      this.panes.prepend(this._identity.root);
      return () => this._identity.root.remove();
    });
    this.effect(() => {
      const source = this.entitySource.value;
      const edit = source?.edit.value;
      if (source === null || source === undefined || edit === undefined || !source.isDraft)
        return;
      this._showIdentity(null, null);
      this.panes.prepend(this._identity.root);
      const changes = edit.onChanged.subscribe(() => void this._derive(source));
      const unguard = source.guard(() => {
        const verdict = this._verdicts.get(source.currentRow.peek()?.id ?? '');
        return verdict?.pending ? 'The structure is still being checked' : verdict?.problem ?? null;
      });
      return () => {
        changes.unsubscribe();
        unguard();
        this._identity.root.remove();
        this._verdicts.clear();
      };
    });
  }

  menuActions(): Action[] {
    // no Bulk edit: it would write past the identity a structure derives
    return [...super.menuActions().filter((action) => action.name !== 'Bulk edit…'),
      {name: 'Open as table', icon: 'table', run: () => void this._openAsTable()}];
  }

  /** The list's rows, structure filters included, in a table view. */
  private async _openAsTable(): Promise<void> {
    const source = this.listSource;
    const query = source.query.peek();
    const parsed = typeof query === 'string' ? Filters.parse(query) : {root: query, problems: []};
    if (parsed.problems.length > 0) {
      grok.shell.warning(parsed.problems[0].message);
      return;
    }
    const search = source.search.peek();
    const sort = source.sort.peek();
    const df = await MolTrackTables.frame('compound', parsed.root, {deleted: source.deleted.peek(),
      ...(search === '' ? {} : {search}), ...(sort === '' ? {} : {sort})});
    grok.shell.addTableView(df);
  }

  private async _derive(source: DomainSource): Promise<void> {
    const row = source.currentRow.peek() as RowView<CompoundRow> | null;
    if (row === null || !Rows.isDraft(row))
      return;
    const input = row.original_structure ?? '';
    if (this._verdicts.get(row.id)?.input === input)
      return;
    const verdict: Verdict = {input, pending: input !== '', problem: null};
    this._verdicts.set(row.id, verdict);
    let identity: CompoundIdentity | null = null;
    if (input !== '') {
      const rule = MolTrackTables.matchRule();
      try {
        identity = await Chemistry.identify(input, rule);
        const existing = await MolTrackTables.registeredAs(identity.hash_mol, rule);
        if (existing !== null) {
          const id = MolTrackTables.corporateId(this.table, 'molregno', existing.molregno);
          verdict.problem = `Already registered as ${id}${existing.common_name ? ` (${existing.common_name})` : ''}`;
        }
      }
      catch (e) {
        verdict.problem = e instanceof Error ? e.message : String(e);
        identity = null;
      }
    }
    if (this._verdicts.get(row.id) !== verdict)
      return;
    verdict.pending = false;
    for (const column of IDENTITY_COLUMNS) {
      const value = identity?.[column] ?? null;
      if (!(isEmpty(value) && isEmpty(row[column])))
        (row as Record<string, unknown>)[column] = value;
    }
    this.form.peek()?.input('original_structure')?.revalidate();
    this._showIdentity(identity, verdict.problem, input);
  }

  private _showIdentity(identity: Partial<CompoundIdentity> | null, problem: string | null,
    original: string | null = null): void {
    const shown = (column: keyof CompoundIdentity) => isEmpty(identity?.[column]) ? '—' : identity![column]!;
    const rule = shown('match_rule');
    const ruleName = `${rule.charAt(0)}${rule.slice(1).toLowerCase().replace(/_/g, ' ')}`;
    const field = (caption: string, value: string | Node) => ObjectForm.readonlyField(caption, caption, value).row;
    const formula = field('Formula', shown('formula'));
    this._identityForm?.dispose();
    const form = this._identityForm = this.runInScope(() => new Form({labelColumn: this.labelColumn}))
      .addElement(formula)
      .addElement(field('InChIKey', shown('inchikey'))).addElement(field('Matching rule', ruleName));
    const structure = identity?.structure;
    // the parent is drawn only where it differs from the Structure the form draws (a salt or a charge taken off)
    if (!isEmpty(structure)) {
      void Chemistry.same(original ?? '', structure!).catch(() => false).then((same) => {
        if (!same)
          formula.before(field('Parent', grok.chem.drawMolecule(structure!, 160, 120)));
      });
    }
    this._identity.body.replaceChildren(...(problem === null ? [] : [ui.divText(problem, 'moltrack-ems-problem')]),
      identity === null ? ui.divText(problem === null ? 'Draw a structure to see its identity' : '',
        'moltrack-ems-hint') : form.root);
  }
}
