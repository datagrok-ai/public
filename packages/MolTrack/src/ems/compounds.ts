import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {DomainFilters} from '@datagrok-libraries/ems/src/dg';
import type {DomainApp, DomainTable} from '@datagrok-libraries/ems/src/dg';
import type {DomainSource} from '@datagrok-libraries/ems';
import {Filters, Form, Rows, Section, isEmpty} from '@datagrok-libraries/u2';
import {ObjectForm} from '@datagrok-libraries/u2/dg';
import type {Action, RowView, Scope} from '@datagrok-libraries/u2';
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

interface IdentityPane {
  app: DomainApp;
  self: Scope;
  section: Section;
  form: Form | null;
  verdicts: Map<string, Verdict>;
}

// EMS has no server write hook: a draft's identity is derived in the browser as its structure is drawn
export class Compounds {
  static readonly DERIVED = ['structure', 'formula', 'inchikey', 'inchi', 'hash_mol', 'hash_all_layers',
    'hash_stereo_insensitive', 'hash_tautomer_insensitive', 'match_rule'];

  /** The Compounds app's `setup`: structure search, the structure's verdict as a validator, the Identity pane. */
  static setup(app: DomainApp, self: Scope): void {
    const pane: IdentityPane = {app, self, section: new Section({title: 'Identity', collapsible: false}), form: null,
      verdicts: new Map()};
    pane.section.root.classList.add('moltrack-identity');
    self.own(DomainFilters.registerOperators(app.table.address, StructureSearch.operators()));
    self.own((app.table as DomainTable<CompoundRow>).validators.add('original_structure', (_, row) =>
      pane.verdicts.get(row.id)?.problem ?? null));
    app.listSource.hosted = (condition, signal) => StructureSearch.matching(condition, signal);
    // a saved compound shows the identity stored with it; only a draft's is derived
    self.effect(() => {
      const source = app.entitySource.value;
      if (source === null || source.isDraft)
        return;
      const row = source.currentRow.value as RowView<CompoundRow> | null;
      if (row === null)
        return;
      Compounds._showIdentity(pane, row, null, row.original_structure);
      app.panes.prepend(pane.section.root);
      return () => Compounds._hideIdentity(pane);
    });
    self.effect(() => {
      const source = app.entitySource.value;
      const edit = source?.edit.value;
      if (source === null || source === undefined || edit === undefined || !source.isDraft)
        return;
      Compounds._showIdentity(pane, null, null);
      app.panes.prepend(pane.section.root);
      const changes = edit.onChanged.subscribe(() => void Compounds._derive(pane, source));
      const unguard = source.guard(() => {
        const verdict = pane.verdicts.get(source.currentRow.peek()?.id ?? '');
        return verdict?.pending ? 'The structure is still being checked' : verdict?.problem ?? null;
      });
      return () => {
        changes.unsubscribe();
        unguard();
        Compounds._hideIdentity(pane);
        pane.verdicts.clear();
      };
    });
  }

  static menu(base: Action[], app: DomainApp): Action[] {
    // no Bulk edit: it would write past the identity a structure derives
    return [...base.filter((action) => action.text !== 'Bulk edit…'),
      {text: 'Open as table', icon: 'table', run: () => void Compounds._openAsTable(app)}];
  }

  /** The list's rows, structure filters included, in a table view. */
  private static async _openAsTable(app: DomainApp): Promise<void> {
    const source = app.listSource;
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

  private static async _derive(pane: IdentityPane, source: DomainSource): Promise<void> {
    const row = source.currentRow.peek() as RowView<CompoundRow> | null;
    if (row === null || !Rows.isDraft(row))
      return;
    const input = row.original_structure ?? '';
    if (pane.verdicts.get(row.id)?.input === input)
      return;
    const verdict: Verdict = {input, pending: input !== '', problem: null};
    pane.verdicts.set(row.id, verdict);
    let identity: CompoundIdentity | null = null;
    if (input !== '') {
      const rule = MolTrackTables.matchRule();
      try {
        identity = await Chemistry.identify(input, rule);
        const existing = await MolTrackTables.registeredAs(identity.hash_mol, rule);
        if (existing !== null) {
          const id = MolTrackTables.corporateId(pane.app.table, 'molregno', existing.molregno);
          verdict.problem = `Already registered as ${id}${existing.common_name ? ` (${existing.common_name})` : ''}`;
        }
      }
      catch (e) {
        verdict.problem = e instanceof Error ? e.message : String(e);
        identity = null;
      }
    }
    if (pane.verdicts.get(row.id) !== verdict)
      return;
    verdict.pending = false;
    for (const column of IDENTITY_COLUMNS) {
      const value = identity?.[column] ?? null;
      if (!(isEmpty(value) && isEmpty(row[column])))
        (row as Record<string, unknown>)[column] = value;
    }
    pane.app.form.peek()?.input('original_structure')?.revalidate();
    Compounds._showIdentity(pane, identity, verdict.problem, input);
  }

  private static _showIdentity(pane: IdentityPane, identity: Partial<CompoundIdentity> | null, problem: string | null,
    original: string | null = null): void {
    const shown = (column: keyof CompoundIdentity) => isEmpty(identity?.[column]) ? '—' : identity![column]!;
    const rule = shown('match_rule');
    const ruleName = `${rule.charAt(0)}${rule.slice(1).toLowerCase().replace(/_/g, ' ')}`;
    const field = (caption: string, value: string | Node) => ObjectForm.readonlyField(caption, caption, value).row;
    const formula = field('Formula', shown('formula'));
    pane.form?.dispose();
    const form = pane.form = pane.self.build(() => new Form({labelColumn: pane.app.labelColumn}))
      .addElement(formula)
      .addElement(field('InChIKey', shown('inchikey'))).addElement(field('Matching rule', ruleName));
    const structure = identity?.structure;
    // the parent is drawn only where it differs from the Structure the form draws (a salt or a charge taken off)
    if (!isEmpty(structure)) {
      form.task(async (t) => {
        const same = await Chemistry.same(original ?? '', structure!).catch(() => false);
        if (!same)
          t.build(() => formula.before(field('Parent', grok.chem.drawMolecule(structure!, 160, 120))));
      });
    }
    pane.section.body.replaceChildren(...(problem === null ? [] : [ui.divText(problem, 'moltrack-problem')]),
      identity === null ? ui.divText(problem === null ? 'Draw a structure to see its identity' : '',
        'moltrack-hint') : form.root);
  }

  private static _hideIdentity(pane: IdentityPane): void {
    pane.section.root.remove();
    pane.form?.dispose();
    pane.form = null;
  }
}
