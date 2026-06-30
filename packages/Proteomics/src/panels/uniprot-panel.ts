import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {fetchUniProtEntry} from '../analysis/subcellular-location';
import {getGroups, GroupAssignment} from '../analysis/experiment-setup';
import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';


/** UniProt REST API response types */
interface UniProtEntry {
  accession?: string;
  proteinDescription?: {
    recommendedName?: { fullName?: { value?: string } };
    submissionNames?: Array<{ fullName?: { value?: string } }>;
  };
  genes?: Array<{ geneName?: { value?: string } }>;
  organism?: { scientificName?: string; commonName?: string };
  comments?: Array<{
    commentType?: string;
    texts?: Array<{ value?: string }>;
  }>;
  uniProtKBCrossReferences?: Array<{
    database?: string;
    id?: string;
    properties?: Array<{ key?: string; value?: string }>;
  }>;
}

interface GoTerms {
  mf: string[];
  bp: string[];
  cc: string[];
}


/**
 * Parse a UniProt accession from a raw protein ID value.
 * Handles:
 *   - Bare accessions: "P04637"
 *   - UniProt format prefixes: "sp|P04637|P53_HUMAN" or "tr|Q67890|..."
 *   - Semicolon-delimited protein groups: "P12345;Q67890" (uses first)
 */
export function parseAccession(raw: string): string | null {
  if (!raw || !raw.trim())
    return null;

  // Take first entry if semicolon-delimited
  const first = raw.split(';')[0].trim();
  if (!first)
    return null;

  // Handle sp|ACC|NAME or tr|ACC|NAME format
  const pipeMatch = first.match(/^(?:sp|tr)\|([A-Za-z0-9_-]+)\|/);
  if (pipeMatch)
    return pipeMatch[1];

  // Handle CPTAC/UPS format: P00167ups|CYB5_HUMAN_UPS — extract accession before "ups" suffix
  const upsMatch = first.match(/^([A-Za-z0-9]+?)ups\|/i);
  if (upsMatch)
    return upsMatch[1];

  // Handle generic ACC|NAME format (single pipe, no sp/tr prefix)
  if (first.includes('|'))
    return first.split('|')[0].trim();

  // Bare accession -- strip trailing "ups" suffix (CPTAC spike-in format) then return
  const bareMatch = first.match(/^([A-Za-z0-9_-]+)$/);
  if (bareMatch) {
    const acc = bareMatch[1];
    return acc.replace(/ups$/i, '');
  }

  return null;
}


/** Fetch protein data from UniProt. Delegates to the shared subcellular-location
 * module's cached fetcher so re-clicking a protein doesn't re-hit the network
 * (closes the folded cache-uniprot todo); the shared fetcher keeps the same
 * fetchProxy / warn+null discipline this panel relied on. */
async function fetchUniProtData(accession: string): Promise<UniProtEntry | null> {
  return (await fetchUniProtEntry(accession)) as UniProtEntry | null;
}


/** Extract and group GO terms by aspect (F/P/C). */
function extractGoTerms(entry: UniProtEntry): GoTerms {
  const mf: string[] = [];
  const bp: string[] = [];
  const cc: string[] = [];

  for (const xref of entry.uniProtKBCrossReferences ?? []) {
    if (xref.database !== 'GO') continue;
    const term = xref.properties?.find((p) => p.key === 'GoTerm')?.value ?? '';
    if (term.startsWith('F:')) mf.push(term.substring(2));
    else if (term.startsWith('P:')) bp.push(term.substring(2));
    else if (term.startsWith('C:')) cc.push(term.substring(2));
  }
  return {mf, bp, cc};
}


/** Get protein name, falling back to submission name if no recommended name. */
function getProteinName(entry: UniProtEntry): string {
  return entry.proteinDescription?.recommendedName?.fullName?.value ??
    entry.proteinDescription?.submissionNames?.[0]?.fullName?.value ??
    'Unknown';
}


/**
 * R3 / D-11 — locate the host DataFrame for a clicked protein in the UniProt panel.
 * Walks grok.shell.tables for tables whose `getGroups()` is set AND whose primary
 * protein-id column resolves the queried accession in some row. Prefers the
 * table behind the current table view when multiple match.
 */
export function findHostDataFrameForProtein(accession: string): DG.DataFrame | null {
  if (!accession) return null;
  const candidates: DG.DataFrame[] = [];
  for (const df of grok.shell.tables) {
    if (!getGroups(df)) continue;
    const idCol = findColumn(df, SEMTYPE.PROTEIN_ID,
      ['primary protein id', 'protein id', 'uniprot', 'accession']);
    if (!idCol) continue;
    for (let i = 0; i < df.rowCount; i++) {
      const raw = idCol.get(i) as string | null;
      if (raw && parseAccession(raw) === accession) {
        candidates.push(df);
        break;
      }
    }
  }
  if (candidates.length === 0) return null;
  const activeDf = grok.shell.tv?.dataFrame ?? null;
  if (activeDf && candidates.includes(activeDf)) return activeDf;
  return candidates[0];
}

/**
 * R3 / D-11 — compact SVG bar chart (one bar per experimental group, mean ± SD
 * whiskers, magenta/cyan per D-04). Returns null when the protein has no row
 * in the host DataFrame; the caller renders the empty-state message instead.
 * Returns the empty-state element when the row exists but every group has zero
 * non-null observations.
 *
 * Source: 14-RESEARCH.md §"Pattern 6 — Inline SVG bar chart in a panel".
 */
export function renderPerGroupBars(df: DG.DataFrame, accession: string): HTMLElement | null {
  const groups = getGroups(df);
  if (!groups) return null;

  const idCol = findColumn(df, SEMTYPE.PROTEIN_ID,
    ['primary protein id', 'protein id', 'uniprot', 'accession']);
  if (!idCol) return null;

  let rowIdx = -1;
  for (let i = 0; i < df.rowCount; i++) {
    const raw = idCol.get(i) as string | null;
    if (raw && parseAccession(raw) === accession) { rowIdx = i; break; }
  }
  if (rowIdx < 0) return null;

  const stats = [groups.group1, groups.group2].map((g) => {
    const vals = g.columns
      .map((n) => df.col(n))
      .filter((c): c is DG.Column => c != null)
      .map((c) => (c.isNone(rowIdx) ? NaN : c.get(rowIdx) as number))
      .filter((v) => typeof v === 'number' && !isNaN(v) && v !== DG.FLOAT_NULL);
    const n = vals.length;
    if (n === 0) return {name: g.name, mean: NaN, sd: NaN, n: 0};
    const mean = vals.reduce((s, v) => s + v, 0) / n;
    const variance = n > 1 ? vals.reduce((s, v) => s + (v - mean) ** 2, 0) / (n - 1) : 0;
    return {name: g.name, mean, sd: Math.sqrt(variance), n};
  });
  if (stats.every((s) => s.n === 0))
    return ui.divText('No per-group quantities available for this protein');

  // SVG layout — locked spacing per 14-UI-SPEC §"Spacing Scale > Per-group bar-chart".
  const W = 200, H = 120, BAR_W = 24, GAP = 8, PAD = 24;
  const maxVal = Math.max(
    ...stats.map((s) => (isNaN(s.mean) ? 0 : s.mean + (isNaN(s.sd) ? 0 : s.sd))),
  );
  const safeMax = maxVal > 0 ? maxVal : 1;
  const scale = (v: number) => H - PAD - (v / safeMax) * (H - 2 * PAD);
  const COLORS = ['#FF00FF', '#00FFFF']; // D-04 LOCKED palette: group1=magenta, group2=cyan
  const svgNs = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNs, 'svg');
  svg.setAttribute('width', String(W));
  svg.setAttribute('height', String(H));

  stats.forEach((s, i) => {
    const x = PAD + i * (BAR_W + GAP * 2);
    if (s.n === 0) return;
    const y = scale(s.mean);
    const rect = document.createElementNS(svgNs, 'rect');
    rect.setAttribute('x', String(x));
    rect.setAttribute('y', String(y));
    rect.setAttribute('width', String(BAR_W));
    rect.setAttribute('height', String(Math.max(0, H - PAD - y)));
    rect.setAttribute('fill', COLORS[i] ?? '#AAAAAA');
    svg.appendChild(rect);

    // Mean±SD whiskers only when SD is defined (n > 1).
    if (s.n > 1 && !isNaN(s.sd)) {
      const whiskerX = x + BAR_W / 2;
      const yTop = scale(s.mean + s.sd);
      const yBot = scale(s.mean - s.sd);
      const line = document.createElementNS(svgNs, 'line');
      line.setAttribute('x1', String(whiskerX));
      line.setAttribute('x2', String(whiskerX));
      line.setAttribute('y1', String(yTop));
      line.setAttribute('y2', String(yBot));
      line.setAttribute('stroke', '#333');
      line.setAttribute('stroke-width', '1');
      svg.appendChild(line);
    }

    const label = document.createElementNS(svgNs, 'text');
    label.setAttribute('x', String(x + BAR_W / 2));
    label.setAttribute('y', String(H - 4));
    label.setAttribute('text-anchor', 'middle');
    label.setAttribute('font-size', '10');
    label.textContent = `${s.name} (n=${s.n})`;
    svg.appendChild(label);
  });

  const container = ui.divV([svg]);
  stats.forEach((s) => {
    if (s.n > 0) {
      const txt = ui.divText(`${s.name}: mean=${s.mean.toFixed(2)}  SD=${s.sd.toFixed(2)}`);
      txt.style.fontSize = '0.85em';
      container.appendChild(txt);
    }
  });
  return container;
}


/** Build the widget DOM for a successfully fetched UniProt entry. */
export function renderUniProtWidget(entry: UniProtEntry, accession: string): HTMLElement {
  const container = document.createElement('div');

  // Prominent link to full UniProt entry
  const link = ui.link(
    `UniProt: ${accession}`,
    `https://www.uniprot.org/uniprot/${accession}`,
    'Open full UniProt entry',
  );
  link.style.fontWeight = 'bold';
  link.style.display = 'block';
  link.style.marginBottom = '8px';
  container.appendChild(link);

  // Summary table
  const proteinName = getProteinName(entry);
  const gene = entry.genes?.[0]?.geneName?.value ?? 'Unknown';
  const organism = entry.organism?.scientificName ?? 'Unknown';

  const funcComment = entry.comments?.find((c) => c.commentType === 'FUNCTION');
  let funcText = funcComment?.texts?.[0]?.value ?? '';
  if (funcText.length > 200)
    funcText = funcText.substring(0, 200) + '...';

  const map: Record<string, string> = {
    'Protein': proteinName,
    'Gene': gene,
    'Organism': organism,
  };
  if (funcText)
    map['Function'] = funcText;

  container.appendChild(ui.tableFromMap(map));

  // GO terms grouped by category
  const go = extractGoTerms(entry);
  const goCategories: Array<{label: string; terms: string[]}> = [
    {label: 'Molecular Function', terms: go.mf},
    {label: 'Biological Process', terms: go.bp},
    {label: 'Cellular Component', terms: go.cc},
  ];

  const hasGoTerms = goCategories.some((c) => c.terms.length > 0);
  if (hasGoTerms) {
    const goHeader = ui.divText('GO Terms');
    goHeader.style.fontWeight = 'bold';
    goHeader.style.marginTop = '8px';
    goHeader.style.marginBottom = '4px';
    container.appendChild(goHeader);

    for (const cat of goCategories) {
      if (cat.terms.length === 0) continue;
      const label = ui.divText(cat.label);
      label.style.fontWeight = 'bold';
      label.style.fontSize = '0.9em';
      label.style.marginTop = '4px';
      container.appendChild(label);

      const termText = cat.terms.slice(0, 5).join(', ');
      const termDiv = ui.divText(termText);
      termDiv.style.fontSize = '0.85em';
      termDiv.style.marginLeft = '8px';
      container.appendChild(termDiv);
    }
  }

  // R3 / D-11 — per-group quantities bar chart (Phase 14)
  // The panel's `accession` parameter is the canonical source; threading it
  // through findHostDataFrameForProtein avoids depending on whichever optional
  // field the UniProtData shape happens to expose.
  const hostDf = findHostDataFrameForProtein(accession);
  if (hostDf) {
    const perGroupHeader = ui.divText('Per-Group Quantities');
    perGroupHeader.style.fontWeight = 'bold';
    perGroupHeader.style.marginTop = '8px';
    perGroupHeader.style.marginBottom = '4px';
    container.appendChild(perGroupHeader);
    const bars = renderPerGroupBars(hostDf, accession);
    container.appendChild(bars ?? ui.divText('No per-group quantities available for this protein'));
  }

  return container;
}


/**
 * Main panel function. Accepts a raw protein ID string, returns a DG.Widget
 * with async loading (spinner while fetching).
 */
export function uniprotPanel(proteinId: string): DG.Widget {
  const accession = parseAccession(proteinId);

  if (!accession)
    return new DG.Widget(ui.divText('No valid UniProt accession found'));

  return new DG.Widget(ui.wait(async () => {
    const data = await fetchUniProtData(accession);
    if (!data)
      return ui.divText(`Unable to fetch UniProt data for ${accession}`);
    return renderUniProtWidget(data, accession);
  }));
}
