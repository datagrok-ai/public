/** Finding a program from the landing page: by its code, its name or an alias, across every
 * portfolio the user can see, and the programs the user opened last. An alias is any short text a
 * registry delivered with the program (its target, or a value in `details` such as a molecule
 * number), so a portfolio needs no configuration to be searchable by it.
 *
 * Pure: the caller reads the programs the user can see and their saved recents. */

export interface FindProgram {
  id: string;
  code: string;
  name: string;
  /** The portfolio's code and name. */
  portfolio: string;
  portfolioName: string;
  status?: string | null;
  target?: string | null;
  details?: Record<string, unknown> | null;
}

export interface FindTerm {
  /** Where the value came from: code, name, target, or the details key. */
  field: string;
  value: string;
}

export interface FindHit {
  program: FindProgram;
  /** The term that matched best. */
  term: FindTerm;
  /** 0 is best: exact, then start of the code, start of a term, start of a word, anywhere. */
  rank: number;
}

/** Fewer characters than this find nothing: one letter matches most of a board. */
export const MIN_QUERY = 2;
export const MAX_HITS = 12;
/** A details value longer than this is prose, not an alias. */
const MAX_ALIAS = 60;

/** Lower case, with dashes, underscores, dots and runs of spaces as one space: "ABC-123", "abc 123"
 * and "abc_123" are the same. */
export const normalize = (s: string): string => s.toLowerCase().replace(/[\s_\-–—·.]+/g, ' ').trim();

/** The terms [p] can be found by. Details values that are numbers or short texts count; nested
 * objects and long texts do not. */
export function findTerms(p: FindProgram): FindTerm[] {
  const terms: FindTerm[] = [{field: 'code', value: p.code}];
  if (p.name && p.name !== p.code)
    terms.push({field: 'name', value: p.name});
  if (p.target)
    terms.push({field: 'target', value: p.target});
  for (const [k, v] of Object.entries(p.details ?? {})) {
    const s = typeof v === 'number' ? String(v) : typeof v === 'string' ? v.trim() : '';
    if (s && s.length <= MAX_ALIAS)
      terms.push({field: k.replace(/_/g, ' '), value: s});
  }
  return terms;
}

/** How well [query] (normalized) matches [value]; null when it does not. */
function rankOf(query: string, value: string, isCode: boolean): number | null {
  const v = normalize(value);
  if (v === query)
    return 0;
  if (v.startsWith(query))
    return isCode ? 1 : 2;
  if (v.split(' ').some((w) => w.startsWith(query)))
    return 3;
  return v.includes(query) ? 4 : null;
}

/** The programs [query] finds, best first, then by portfolio and code; at most [limit]. */
export function findPrograms(query: string, programs: FindProgram[], limit = MAX_HITS): FindHit[] {
  const q = normalize(query);
  if (q.length < MIN_QUERY)
    return [];
  const hits: FindHit[] = [];
  for (const p of programs) {
    let best: FindHit | null = null;
    for (const term of findTerms(p)) {
      const rank = rankOf(q, term.value, term.field === 'code');
      if (rank !== null && (!best || rank < best.rank))
        best = {program: p, term, rank};
    }
    if (best)
      hits.push(best);
  }
  return hits.sort((a, b) => a.rank - b.rank || a.program.portfolio.localeCompare(b.program.portfolio) ||
    a.program.code.localeCompare(b.program.code)).slice(0, limit);
}

/** Why [hit] was found, for the result line: nothing when the code or name says it. */
export function matchedBy(hit: FindHit): string | null {
  return hit.term.field === 'code' || hit.term.field === 'name' ? null : `${hit.term.field}: ${hit.term.value}`;
}

/** One program the user opened, by id (a code may change), newest first in the list. */
export interface RecentProgram {
  id: string;
  portfolio: string;
  at: string;
}

export const MAX_RECENT = 8;

/** The saved recents; anything unreadable is dropped, never thrown. */
export function parseRecents(json: string | null | undefined): RecentProgram[] {
  try {
    const list = JSON.parse(json ?? '[]');
    return Array.isArray(list) ? list.filter((r) => r && typeof r.id === 'string' && typeof r.portfolio === 'string')
      .map((r) => ({id: r.id, portfolio: r.portfolio, at: String(r.at ?? '')})) : [];
  } catch {
    return [];
  }
}

/** [list] with [entry] first and no other entry for the same program, at most [max] long. */
export function addRecent(list: RecentProgram[], entry: RecentProgram, max = MAX_RECENT): RecentProgram[] {
  return [entry, ...list.filter((r) => r.id !== entry.id)].slice(0, max);
}

/** The recents the user can still see, in order: a program gone or no longer shared is skipped. */
export function resolveRecents(recents: RecentProgram[], programs: FindProgram[]): FindProgram[] {
  const byId = new Map(programs.map((p) => [p.id, p]));
  return recents.map((r) => byId.get(r.id)).filter((p): p is FindProgram => !!p);
}
