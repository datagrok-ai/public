import * as DG from 'datagrok-api/dg';
import type {DomainImporter, DomainImportVerdict, DomainTable} from '@datagrok-libraries/ems/src/dg';
import {isEmpty} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, HASH_COLUMNS, IDENTITY_COLUMNS} from './chemistry';
import type {CompoundIdentity, MatchRule} from './chemistry';
import {MolTrackTables} from './tables';

const LOOKUP_CHUNK = 500;

// bulk registration is the stock import of Compounds: this is the chemistry it runs on every row first
export class BulkRegistration {
  static importer(compounds: DomainTable<CompoundRow>): DomainImporter {
    return {derives: IDENTITY_COLUMNS,
      prepare: (rows, progress) => BulkRegistration.prepare(rows, compounds, progress)};
  }

  static async prepare(rows: Record<string, unknown>[], compounds: DomainTable<CompoundRow>,
    progress: (done: number) => void = () => {}): Promise<DomainImportVerdict[]> {
    const rule = MolTrackTables.matchRule();
    const identities: (CompoundIdentity | null)[] = [];
    for (const [i, row] of rows.entries()) {
      const input = isEmpty(row.original_structure) ? '' : String(row.original_structure);
      identities.push(input === '' ? null : await Chemistry.identify(input, rule).catch(() => null));
      progress(i + 1);
    }
    const registered = await BulkRegistration._registered(
      [...new Set(identities.flatMap((identity) => identity === null ? [] : [identity.hash_mol]))], rule);
    const first = new Map<string, number>();
    return identities.map((identity, i): DomainImportVerdict => {
      if (identity === null)
        return {error: isEmpty(rows[i].original_structure) ? 'No structure' : 'Cannot read the structure'};
      const existing = registered.get(identity.hash_mol);
      if (existing !== undefined) {
        return {existing: existing.id,
          note: `Already registered as ${MolTrackTables.corporateId(compounds, 'molregno', existing.molregno)}`};
      }
      const seen = first.get(identity.hash_mol);
      if (seen !== undefined)
        return {sameAs: seen, note: `Same compound as row ${seen + 1}`};
      first.set(identity.hash_mol, i);
      return {values: Object.fromEntries(IDENTITY_COLUMNS.map((column) => [column, identity[column]]))};
    });
  }

  // under the matching rule in force: a compound registered under another rule is not "the same" here
  private static async _registered(hashes: string[], rule: MatchRule): Promise<Map<string, CompoundRow>> {
    const column = HASH_COLUMNS[rule];
    const found = new Map<string, CompoundRow>();
    for (let at = 0; at < hashes.length; at += LOOKUP_CHUNK) {
      const filter = DG.cond(column, '=', hashes.slice(at, at + LOOKUP_CHUNK));
      const rows = await moltrackEmsDb.compounds.query({filter, columns: [column, 'molregno'], sort: 'molregno'});
      for (const row of rows) {
        if (!found.has(row[column]!))
          found.set(row[column]!, row);
      }
    }
    return found;
  }
}
