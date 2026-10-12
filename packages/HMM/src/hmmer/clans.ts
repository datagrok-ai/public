// Synced from Rusty-HMMER web/clans.ts by web/sync-datagrok.mjs — do not edit here.
/** Pfam clan competition: overlapping matches of families in the same clan
 * are resolved in favour of the most significant one, as Pfam reports them.
 *
 * The rule was established black-box against Pfam's own pfam_scan.pl (run
 * unmodified as an oracle, web/pfam/clan-oracle.sh; it is not translated) and
 * is checked against it by web/pfam/clans.test.ts:
 * - the domains of one sequence are ranked by their E-value as hmmscan prints
 *   it (`%.2g`, so ties are common), ties by the rank of their hit in
 *   hmmscan's hit list, then by domain order;
 * - in rank order, a domain is dropped when its alignment (ali from..to,
 *   inclusive) overlaps a domain already kept whose family is in the same clan
 *   (the same family included); families without a clan are never dropped.
 * Pfam lets some families nest inside others (Pfam-A.hmm.dat `NE` lines);
 * no such pair is in the same clan in the built-in library, and nesting is
 * not modelled for other models. */

export interface ClanDomain {
  /** Pfam clan accession of the domain's family, or null. */
  clan: string | null;
  aliFrom: number;
  aliTo: number;
  /** Independent E-value as hmmscan prints it (`DomainHit.printed.iEvalue`). */
  evalue: number;
  /** Rank of the domain's hit in hmmscan's (sorted) hit list. */
  hitRank: number;
  /** Index of the domain within its hit. */
  domain: number;
}

/** The domains of one sequence that survive clan competition, in input order. */
export function resolveClanOverlaps<T extends ClanDomain>(domains: readonly T[]): T[] {
  const ranked = domains.map((d, i) => ({d, i})).sort((a, b) =>
    a.d.evalue - b.d.evalue || a.d.hitRank - b.d.hitRank || a.d.domain - b.d.domain || a.i - b.i);
  const kept: {d: T; i: number}[] = [];
  for (const entry of ranked) {
    const {clan, aliFrom, aliTo} = entry.d;
    const loses = clan !== null &&
      kept.some(({d}) => d.clan === clan && aliFrom <= d.aliTo && d.aliFrom <= aliTo);
    if (!loses) kept.push(entry);
  }
  return kept.sort((a, b) => a.i - b.i).map(({d}) => d);
}
