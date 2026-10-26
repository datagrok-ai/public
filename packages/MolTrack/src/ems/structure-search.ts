import * as DG from 'datagrok-api/dg';
import {Filters, KIND} from '@datagrok-libraries/u2';
import type {FilterCondition, FilterGroup} from '@datagrok-libraries/u2';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, HASH_COLUMNS} from './chemistry';
import {MolTrackTables} from './tables';

// answered in the browser with Crux: the server has no chemistry
export class StructureSearch {
  static readonly SUBSTRUCTURE = 'has substructure';
  static readonly SUPERSTRUCTURE = 'is substructure of';
  static readonly SIMILAR = 'is similar to';
  static readonly EXACT = 'is exactly';
  static readonly PAGE = 5000;

  static readonly OPERATORS = [StructureSearch.SUBSTRUCTURE, StructureSearch.SUPERSTRUCTURE, StructureSearch.SIMILAR,
    StructureSearch.EXACT];

  static register(): () => void {
    return Filters.operators.register(StructureSearch.OPERATORS.map((id) => ({id, label: id, arity: 1 as const,
      kinds: [KIND.STRING], semType: DG.SEMTYPE.MOLECULE, exclusive: true, editor: 'default' as const, host: true})));
  }

  static threshold(): number {
    const value = Number(Chemistry.package().settings?.['SimilarityThreshold']);
    return Number.isFinite(value) && value > 0 && value <= 1 ? value : 0.5;
  }

  /** The query with every structure condition replaced by the ids of the compounds it matches. */
  static async resolve(root: FilterGroup): Promise<FilterGroup> {
    return {...root, nodes: await Promise.all(root.nodes.map(async (n) => Filters.isGroup(n) ?
      StructureSearch.resolve(n) : StructureSearch.OPERATORS.includes(n.operator) ?
        Filters.cond('id', 'in', await StructureSearch.matching(n)) : n))};
  }

  static async matching(condition: FilterCondition, signal?: AbortSignal): Promise<string[]> {
    const query = String(condition.value);
    if (condition.operator === StructureSearch.EXACT) {
      const rule = MolTrackTables.matchRule();
      const {hash_mol: hash} = await Chemistry.identify(query, rule);
      const rows = await moltrackEmsDb.compounds.query({filter: DG.cond(HASH_COLUMNS[rule], '=', hash),
        columns: [HASH_COLUMNS[rule]], deleted: 'include'});
      return rows.map((r) => r.id);
    }
    const {ids, structures} = await StructureSearch._column(condition.property, signal);
    if (condition.operator === StructureSearch.SIMILAR) {
      const scores = await Chemistry.similarity(structures, query);
      const threshold = StructureSearch.threshold();
      return ids.filter((_, i) => scores[i] >= threshold);
    }
    const hits = condition.operator === StructureSearch.SUBSTRUCTURE ?
      await Chemistry.substructure(structures, query) : await Chemistry.superstructure(structures, query);
    return hits.map((i) => ids[i]);
  }

  private static async _column(property: string, signal?: AbortSignal):
    Promise<{ids: string[], structures: string[]}> {
    const column = property as 'structure' | 'original_structure';
    const ids: string[] = [];
    const structures: string[] = [];
    for (let offset = 0; ; offset += StructureSearch.PAGE) {
      signal?.throwIfAborted();
      const rows = await moltrackEmsDb.compounds.query({columns: [column], sort: 'molregno',
        limit: StructureSearch.PAGE, offset, deleted: 'include'});
      for (const row of rows) {
        ids.push(row.id);
        structures.push(row[column] ?? '');
      }
      if (rows.length < StructureSearch.PAGE)
        return {ids, structures};
    }
  }
}
