import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {DomainApp, DomainFilters} from '@datagrok-libraries/ems/src/dg';
import type {DomainAppOptions, DomainTable} from '@datagrok-libraries/ems/src/dg';
import type {DomainSource} from '@datagrok-libraries/ems';
import {Filters, Rows, Section, isEmpty} from '@datagrok-libraries/u2';
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
  private readonly _identity = ui.divV([], 'moltrack-ems-identity');
  private readonly _identityHeading =
    this.runInScope(() => new Section({title: 'Identity', collapsible: false})).header;

  constructor(options: DomainAppOptions) {
    super(options);
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
      this.panes.prepend(this._identity);
      return () => this._identity.remove();
    });
    this.effect(() => {
      const source = this.entitySource.value;
      const edit = source?.edit.value;
      if (source === null || source === undefined || edit === undefined || !source.isDraft)
        return;
      this._showIdentity(null, null);
      this.panes.prepend(this._identity);
      const changes = edit.onChanged.subscribe(() => void this._derive(source));
      const unguard = source.guard(() => {
        const verdict = this._verdicts.get(source.currentRow.peek()?.id ?? '');
        return verdict?.pending ? 'The structure is still being checked' : verdict?.problem ?? null;
      });
      return () => {
        changes.unsubscribe();
        unguard();
        this._identity.remove();
        this._verdicts.clear();
      };
    });
    // the identity's captions take the form's label column, which the form sets on its own root
    this.effect(() => {
      const form = this.form.value;
      if (form === null)
        return;
      const align = () => this._identity.style.setProperty('--u2-form-label-width',
        form.root.querySelector<HTMLElement>('.u2-form')?.style.getPropertyValue('--u2-form-label-width') ?? '');
      align();
      const observer = new MutationObserver(align);
      observer.observe(form.root, {subtree: true, childList: true, attributes: true, attributeFilter: ['style']});
      return () => observer.disconnect();
    });
  }

  menuActions(): Action[] {
    return [...super.menuActions(), {name: 'Open as table', icon: 'table', run: () => void this._openAsTable()}];
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
    // the parent is drawn, not spelled out — only where it is not the Structure the form draws
    // already (a salt or a charge taken off)
    const parent = ui.div([]);
    const rows = ui.tableFromMap({'Parent': parent, 'Formula': shown('formula'), 'InChIKey': shown('inchikey'),
      'Matching rule': `${rule.charAt(0)}${rule.slice(1).toLowerCase().replace(/_/g, ' ')}`});
    const parentRow = parent.closest('tr');
    if (parentRow !== null)
      parentRow.hidden = true;
    const structure = identity?.structure;
    if (!isEmpty(structure)) {
      void Chemistry.same(original ?? '', structure!).catch(() => false).then((same) => {
        if (same)
          return;
        parent.append(grok.chem.drawMolecule(structure!, 160, 120));
        if (parentRow !== null)
          parentRow.hidden = false;
      });
    }
    this._identity.replaceChildren(this._identityHeading,
      ...(problem === null ? [] : [ui.divText(problem, 'moltrack-ems-problem')]),
      identity === null ? ui.divText(problem === null ? 'Draw a structure to see its identity' : '',
        'moltrack-ems-hint') : rows);
  }
}
