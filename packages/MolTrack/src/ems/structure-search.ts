import * as DG from 'datagrok-api/dg';
import {DomainQueries} from '@datagrok-libraries/ems';
import type {DomainTable} from '@datagrok-libraries/ems/src/dg';
import type {FilterCondition} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, HASH_COLUMNS} from './chemistry';
import {MolTrackTables} from './tables';

/** The filter builder's structure conditions, answered by Crux in the browser: the server has no chemistry. */
export class StructureSearch {
  static readonly SUBSTRUCTURE = 'has substructure';
  static readonly SUPERSTRUCTURE = 'is substructure of';
  static readonly SIMILAR = 'is similar to';
  static readonly EXACT = 'is exactly';
  static readonly THRESHOLD = 0.7;
  static readonly PAGE = 5000;

  static register(compounds: DomainTable<CompoundRow>): void {
    for (const id of [StructureSearch.SUBSTRUCTURE, StructureSearch.SUPERSTRUCTURE, StructureSearch.SIMILAR,
      StructureSearch.EXACT]) {
      compounds.operators.add({id, label: id, semType: DG.SEMTYPE.MOLECULE,
        rows: (c, signal) => StructureSearch.matching(c, signal),
        ...(id === StructureSearch.SIMILAR ? {parameter: {name: 'threshold', label: 'Similarity at least',
          default: StructureSearch.THRESHOLD, min: 0, max: 1, step: 0.05}} : {})});
    }
  }

  static async matching(condition: FilterCondition, signal?: AbortSignal): Promise<string[]> {
    const query = String(condition.value);
    if (condition.operator === StructureSearch.EXACT) {
      const rule = MolTrackTables.matchRule();
      const {hash_mol: hash} = await Chemistry.identify(query, rule);
      const rows = await moltrackEmsDb.compounds.query({filter: DG.cond(HASH_COLUMNS[rule], '=', hash),
        columns: [HASH_COLUMNS[rule]]});
      return rows.map((r) => r.id);
    }
    const {ids, structures} =
      await StructureSearch._column(condition.property as 'structure' | 'original_structure', signal);
    signal?.throwIfAborted();
    if (condition.operator === StructureSearch.SIMILAR) {
      const threshold = condition.options?.threshold;
      const scores = await Chemistry.similarity(structures, query);
      return ids.filter((_, i) => scores[i] >= (typeof threshold === 'number' ? threshold : StructureSearch.THRESHOLD));
    }
    const hits = condition.operator === StructureSearch.SUBSTRUCTURE ?
      await Chemistry.substructure(structures, query) : await Chemistry.superstructure(structures, query);
    return hits.map((i) => ids[i]);
  }

  private static async _column(column: 'structure' | 'original_structure', signal?: AbortSignal):
    Promise<{ids: string[], structures: string[]}> {
    const rows = await DomainQueries.all((page) => {
      signal?.throwIfAborted();
      return moltrackEmsDb.compounds.query({columns: [column], sort: 'molregno', ...page});
    }, StructureSearch.PAGE);
    return {ids: rows.map((r) => r.id), structures: rows.map((r) => r[column] ?? '')};
  }
}
