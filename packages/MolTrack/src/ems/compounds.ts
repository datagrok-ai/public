import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {DomainApp} from '@datagrok-libraries/ems/src/dg';
import type {DomainAppOptions, DomainTable} from '@datagrok-libraries/ems/src/dg';
import type {DomainSource} from '@datagrok-libraries/ems';
import {Rows, isEmpty} from '@datagrok-libraries/u2';
import type {Action, RowView} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {Chemistry, IDENTITY_COLUMNS} from './chemistry';
import type {CompoundIdentity} from './chemistry';
import {BulkRegistration} from './registration';
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

  private static readonly _verdicts = new Map<string, Verdict>();
  private readonly _drafts = new Set<string>();
  private readonly _identity = ui.divV([], 'moltrack-ems-identity');

  static configure(compounds: DomainTable<CompoundRow>): void {
    compounds.validators.add('original_structure', (_, row) => CompoundsApp._verdicts.get(row.id)?.problem ?? null);
  }

  constructor(options: DomainAppOptions) {
    super(options);
    this.own(StructureSearch.register());
    this.listSource.hosted = (condition, signal) => StructureSearch.matching(condition, signal);
    this.effect(() => {
      const source = this.entitySource.value;
      const edit = source?.edit.value;
      if (source === null || source === undefined || edit === undefined || !source.isDraft)
        return;
      this._showIdentity(null, null);
      this.panes.prepend(this._identity);
      const changes = edit.onChanged.subscribe(() => void this._derive(source));
      const unguard = source.guard(() => {
        const verdict = CompoundsApp._verdicts.get(source.currentRow.peek()?.id ?? '');
        return verdict?.pending ? 'The structure is still being checked' : verdict?.problem ?? null;
      });
      return () => {
        changes.unsubscribe();
        unguard();
        this._identity.remove();
        for (const id of this._drafts)
          CompoundsApp._verdicts.delete(id);
        this._drafts.clear();
      };
    });
  }

  // no Bulk edit: it would write past the identity a structure derives
  menuActions(): Action[] {
    return super.menuActions().filter((action) => action.name !== 'Bulk edit…')
      .map((action) => action.name !== 'Import…' ? action : {...action,
        run: () => void BulkRegistration.run(this.table as DomainTable<CompoundRow>)
          .then((landed) => landed ? this.listSource.refresh() : null)});
  }

  private async _derive(source: DomainSource): Promise<void> {
    const row = source.currentRow.peek() as RowView<CompoundRow> | null;
    if (row === null || !Rows.isDraft(row))
      return;
    const input = row.original_structure ?? '';
    if (CompoundsApp._verdicts.get(row.id)?.input === input)
      return;
    const verdict: Verdict = {input, pending: input !== '', problem: null};
    CompoundsApp._verdicts.set(row.id, verdict);
    this._drafts.add(row.id);
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
    if (CompoundsApp._verdicts.get(row.id) !== verdict)
      return;
    verdict.pending = false;
    for (const column of IDENTITY_COLUMNS) {
      const value = identity?.[column] ?? null;
      if (!(isEmpty(value) && isEmpty(row[column])))
        (row as Record<string, unknown>)[column] = value;
    }
    this.form.peek()?.input('original_structure')?.revalidate();
    this._showIdentity(identity, verdict.problem);
  }

  private _showIdentity(identity: CompoundIdentity | null, problem: string | null): void {
    this._identity.replaceChildren(ui.h3('Identity'),
      ...(problem === null ? [] : [ui.divText(problem, 'moltrack-ems-problem')]),
      identity === null ? ui.divText(problem === null ? 'Draw a structure to see its identity' : '') :
        ui.divH([grok.chem.drawMolecule(identity.structure, 160, 120), ui.tableFromMap({'Parent': identity.structure,
          'Formula': identity.formula, 'InChIKey': identity.inchikey, 'Matching rule': identity.match_rule})]));
  }
}
