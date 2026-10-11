/** Align to profile: `hmmalign` of a sequence column to one profile HMM, as a
 * Bio MSA engine (`meta.role: sequenceMSA`). Every column of the result is a
 * model position (consensus columns, numbered by model node) or an insertion
 * between two positions, so alignments to the same model are comparable
 * across tables. Traces are computed in parallel workers and the alignment is
 * assembled once (hmmer/tools.ts), identical to C hmmalign. */
import * as DG from 'datagrok-api/dg';
import {HitFlags} from './hmmer/engine.ts';
import type {ProfileAlignment} from './hmmer/tools.ts';
import {HmmerPool} from './pool';
import {columnResidues} from './sequences';
import {fetchPfam, loadBuiltIn, loadFile, loadGermlines, type LoadedModels} from './models';
import {PROFILE_TAG} from './logo';

export const PROFILES = ['Best Pfam family (auto)', 'Pfam family or accession', 'Antibody/TCR germline (auto)',
  'HMM file'];

/** The model most often scoring best over `sequences` (ties: higher total score). */
async function bestModel(pool: HmmerPool, models: LoadedModels, sequences: {name: string; residues: string}[],
  what: string): Promise<number> {
  const results = await pool.scan(models.key, sequences.filter((s) => s.residues));
  const votes = new Map<number, {count: number; score: number}>();
  for (const result of results) {
    const best = result.hits.find((h) => h.flags & HitFlags.reported);
    if (!best) continue;
    const vote = votes.get(best.target) ?? {count: 0, score: 0};
    vote.count++;
    vote.score += best.score;
    votes.set(best.target, vote);
  }
  const ranked = [...votes.entries()].sort(([, a], [, b]) => b.count - a.count || b.score - a.score);
  if (ranked.length === 0) throw new Error(`No ${what} matches these sequences`);
  return ranked[0][0];
}

/** The profile chosen by the engine's parameters: model data and index. */
export async function chooseProfile(pool: HmmerPool, pkg: DG.Package, profile: string, family: string | null,
  hmmFile: DG.FileInfo | null, sequences: {name: string; residues: string}[]):
  Promise<{models: LoadedModels; index: number}> {
  switch (profile) {
  case PROFILES[1]: {
    const wanted = (family ?? '').trim();
    if (!wanted) throw new Error('Enter a Pfam family name (V-set) or accession (PF07686)');
    const library = await loadBuiltIn(pool, pkg);
    const lower = wanted.toLowerCase();
    const index = library.info.models.findIndex((m) =>
      m.name.toLowerCase() === lower || m.accession.split('.')[0].toLowerCase() === lower);
    if (index >= 0) return {models: library, index};
    const accession = /^PF\d{5}$/i.exec(wanted)?.[0];
    if (!accession) throw new Error(`${wanted} is not in the built-in library; give its Pfam accession (PF00000)`);
    return {models: await fetchPfam(pool, [accession.toUpperCase()]), index: 0};
  }
  case PROFILES[2]: {
    const germlines = await loadGermlines(pool, pkg);
    return {models: germlines, index: await bestModel(pool, germlines, sequences, 'antibody or TCR germline model')};
  }
  case PROFILES[3]: {
    if (!hmmFile) throw new Error('Choose an HMM file');
    const models = await loadFile(pool, hmmFile);
    if (models.info.models.length !== 1)
      throw new Error(`${hmmFile.name} has ${models.info.models.length} models; hmmalign aligns to exactly one`);
    return {models, index: 0};
  }
  default: {
    const library = await loadBuiltIn(pool, pkg);
    return {models: library, index: await bestModel(pool, library, sequences, 'Pfam family of the built-in library')};
  }
  }
}

/** Position names of the alignment columns: model node numbers, and `j.n` for
 * the n-th insertion column after node j. */
export function positionNames(consensus: boolean[]): string[] {
  const names: string[] = [];
  let node = 0;
  let insert = 0;
  for (const c of consensus) {
    if (c) {
      names.push(String(++node));
      insert = 0;
    } else
      names.push(`${node}.${++insert}`);
  }
  return names;
}

/** hmmalign of `column` to model `index` of `models`; the aligned Bio column. */
export async function alignColumn(pool: HmmerPool, column: DG.Column<string>, models: LoadedModels, index: number,
  trim: boolean, sequences?: {name: string; residues: string}[]): Promise<DG.Column<string>> {
  sequences ??= (await columnResidues(column)).residues.map((residues, i) => ({name: String(i), residues}));
  if (sequences.length === 0) throw new Error('No sequences to align');
  const model = models.info.models.length === 1 ? models.bytes : await pool.extract(models.bytes, index);
  const aligned: ProfileAlignment = await pool.align(model, 0, sequences, trim);
  const result = DG.Column.fromStrings(`${column.name} (hmmalign)`,
    aligned.rows.map((row) => row.toUpperCase().replace(/\./g, '-')));
  const info = models.info.models[index];
  result.semType = DG.SEMTYPE.MACROMOLECULE;
  result.meta.units = 'fasta';
  result.setTag('aligned', 'SEQ.MSA');
  result.setTag('alphabet', models.info.alphabet === 'amino' ? 'PT' : models.info.alphabet === 'dna' ? 'DNA' : 'RNA');
  result.setTag('.positionNames', positionNames(aligned.consensus).join(', '));
  result.setTag(PROFILE_TAG, [info.name, info.accession].filter((s) => s).join(' '));
  result.setTag(DG.Tags.CellRenderer, 'sequence');
  return result;
}
