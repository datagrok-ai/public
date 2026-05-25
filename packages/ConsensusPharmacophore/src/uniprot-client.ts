/*
 * Thin UniProtKB REST client for the target-lookup feature.
 *
 * The consensus pharmacophore pipeline runs on PDB IDs, but users typically
 * think in terms of "give me the EGFR consensus" — not "give me the consensus
 * for [1XKK, 3W2S, 4WKQ, 5HG5, 5HG8]". This module turns a UniProt accession
 * (e.g. P00533) OR a gene/protein name (e.g. "egfr") into the list of PDB
 * cross-references for that entry. When the user types a name that maps to
 * multiple species (human EGFR, mouse Egfr, fly Egfr, ...), the orchestrator
 * surfaces a picker dialog so they choose one accession before resolving.
 *
 * APIs used:
 *   - GET  https://rest.uniprot.org/uniprotkb/search?query=...&format=json&fields=...
 *   - GET  https://rest.uniprot.org/uniprotkb/{accession}?format=json&fields=...
 *
 * All HTTP goes through `grok.dapi.fetchProxy` so CORS is handled by the
 * Datagrok server and the rest of the package doesn't need a CORS-friendly
 * proxy of its own.
 */

import * as grok from 'datagrok-api/grok';


/** One row of the UniProtKB search results. */
export interface UniProtHit {
  /** Primary UniProt accession (e.g. "P00533"). */
  accession: string;
  /** Organism display name (e.g. "Homo sapiens"). */
  organism: string;
  /** NCBI taxonomy ID (e.g. 9606 for human). */
  taxonomyId: number;
  /** Top gene name (e.g. "EGFR"). May be empty for entries without one. */
  geneName: string;
  /** Recommended protein full name (e.g. "Epidermal growth factor receptor"). */
  proteinName: string;
  /** True for SwissProt (manually reviewed) entries; false for TrEMBL. */
  reviewed: boolean;
  /** Best-effort count of PDB cross-references from the search response.
   *  Useful to show in the picker so the user can prefer the species with
   *  the most structures. May be 0 if the search response didn't expose it
   *  (we then resolve it lazily on accession pick). */
  pdbCount: number;
}


// UniProt accession regex (UniProtKB / SwissProt). Matches the two valid
// shapes per https://www.uniprot.org/help/accession_numbers:
//   [OPQ][0-9][A-Z0-9]{3}[0-9]                — older accessions, e.g. P00533
//   [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2} — newer accessions, e.g. A0A024R161
const ACCESSION_RE =
  /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$/i;

/** True if the query string is a syntactically valid UniProt accession.
 *  Doesn't guarantee the accession exists in UniProtKB — that's checked at
 *  fetch time. Used by the orchestrator to skip the search step when the
 *  user pastes a known accession. */
export function isUniProtAccession(query: string): boolean {
  return ACCESSION_RE.test(query.trim());
}


/** Pull the canonical fields from one element of `json.results`. The UniProt
 *  search API's response shape is verbose; this normalises it into the small
 *  `UniProtHit` shape the rest of the package uses.
 *
 *  Defensive everywhere — `?.` / `??` because UniProt occasionally returns
 *  entries with missing organism / gene-names blocks. */
function _parseUniProtSearchResult(raw: any): UniProtHit | null {
  const accession = raw?.primaryAccession;
  if (!accession) return null;
  const organism = raw?.organism?.scientificName ?? '';
  const taxonomyId = Number(raw?.organism?.taxonId ?? 0) || 0;
  const geneName = (raw?.genes?.[0]?.geneName?.value as string | undefined) ?? '';
  const proteinName = (raw?.proteinDescription?.recommendedName?.fullName?.value as string | undefined) ??
                      (raw?.proteinDescription?.submissionNames?.[0]?.fullName?.value as string | undefined) ??
                      '';
  const reviewed = String(raw?.entryType ?? '').toLowerCase().includes('swiss-prot');
  const pdbCount = ((raw?.uniProtKBCrossReferences ?? []) as any[])
    .filter((r) => r?.database === 'PDB').length;
  return {accession, organism, taxonomyId, geneName, proteinName, reviewed, pdbCount};
}


/**
 * Search UniProtKB for the query string. Returns ranked hits. SwissProt
 * (reviewed) entries are preferred — we add `+AND+reviewed:true` so users get
 * the curated entries first; if they want unreviewed TrEMBL hits they can
 * paste the accession directly.
 *
 *   `limit` — UniProt's `size` parameter, capped at 50 server-side.
 */
export async function searchUniProt(query: string, limit = 25): Promise<UniProtHit[]> {
  const q = query.trim();
  if (!q) return [];
  // Bias the search toward the gene_name and protein_name fields so a query
  // like "egfr" doesn't pull every entry that mentions EGFR in a CC line.
  // The `(field:val)` syntax is documented at
  // https://www.uniprot.org/help/text-search.
  const queryStr = `(gene:${q} OR protein_name:${q}) AND reviewed:true`;
  const fields = [
    'accession', 'organism_name', 'organism_id', 'gene_names',
    'protein_name', 'reviewed', 'xref_pdb',
  ].join(',');
  const url = `https://rest.uniprot.org/uniprotkb/search?` +
    `query=${encodeURIComponent(queryStr)}&format=json&size=${Math.min(limit, 50)}&` +
    `fields=${encodeURIComponent(fields)}`;

  let resp: Response;
  try {
    resp = await grok.dapi.fetchProxy(url);
  } catch (e: any) {
    throw new Error(`UniProt search failed (network): ${e?.message ?? e}`);
  }
  if (!resp.ok) {
    const text = await resp.text().catch(() => '');
    throw new Error(`UniProt search HTTP ${resp.status}: ${text.slice(0, 200)}`);
  }
  let json: any;
  try {
    json = await resp.json();
  } catch (e: any) {
    throw new Error(`UniProt search returned non-JSON: ${e?.message ?? e}`);
  }
  const results: any[] = json?.results ?? [];
  return results
    .map(_parseUniProtSearchResult)
    .filter((h): h is UniProtHit => h !== null);
}


/**
 * Fetch a single UniProt entry by accession. Returns the same `UniProtHit`
 * shape used by `searchUniProt` so the picker UI doesn't have to special-case
 * "user pasted an accession directly".
 */
export async function fetchUniProtById(accession: string): Promise<UniProtHit | null> {
  const acc = accession.trim().toUpperCase();
  if (!isUniProtAccession(acc))
    throw new Error(`Not a UniProt accession: ${accession}`);
  const fields = [
    'accession', 'organism_name', 'organism_id', 'gene_names',
    'protein_name', 'reviewed', 'xref_pdb',
  ].join(',');
  const url = `https://rest.uniprot.org/uniprotkb/${acc}?format=json&` +
    `fields=${encodeURIComponent(fields)}`;

  let resp: Response;
  try {
    resp = await grok.dapi.fetchProxy(url);
  } catch (e: any) {
    throw new Error(`UniProt accession fetch failed (network): ${e?.message ?? e}`);
  }
  if (resp.status === 404) return null;
  if (!resp.ok) {
    const text = await resp.text().catch(() => '');
    throw new Error(`UniProt accession fetch HTTP ${resp.status}: ${text.slice(0, 200)}`);
  }
  const json = await resp.json();
  return _parseUniProtSearchResult(json);
}


/**
 * Return the (deduplicated, alphabetically sorted) PDB IDs cross-referenced
 * from a UniProt entry. Uses the entry-fetch endpoint with `fields=xref_pdb`
 * so we don't pull the full entry payload (which can be ~MB for well-studied
 * proteins like EGFR).
 */
export async function getPdbIdsForAccession(accession: string): Promise<string[]> {
  const acc = accession.trim().toUpperCase();
  if (!isUniProtAccession(acc))
    throw new Error(`Not a UniProt accession: ${accession}`);
  const url = `https://rest.uniprot.org/uniprotkb/${acc}?format=json&fields=xref_pdb`;

  let resp: Response;
  try {
    resp = await grok.dapi.fetchProxy(url);
  } catch (e: any) {
    throw new Error(`UniProt PDB-xref fetch failed (network): ${e?.message ?? e}`);
  }
  if (resp.status === 404) return [];
  if (!resp.ok) {
    const text = await resp.text().catch(() => '');
    throw new Error(`UniProt PDB-xref fetch HTTP ${resp.status}: ${text.slice(0, 200)}`);
  }
  const json = await resp.json();
  const refs: string[] = ((json?.uniProtKBCrossReferences ?? []) as any[])
    .filter((r) => r?.database === 'PDB')
    .map((r) => String(r?.id ?? '').toUpperCase())
    .filter((s) => /^[1-9][A-Z0-9]{3}$/.test(s));
  return [...new Set(refs)].sort();
}
