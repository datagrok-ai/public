/** Pure PDB file header parser. All fields are optional to handle
 *  files from RCSB, MOE, Schrödinger, and other sources. */

export interface PdbEntity {
  molId: string;
  molecule?: string;
  chains?: string[];
}

export interface PdbLigand {
  id: string;
  chain?: string;
  name?: string;
  formula?: string;
}

export interface PdbCitation {
  title?: string;
  journal?: string;
  doi?: string;
  pmid?: number;
}

export interface PdbHeaderInfo {
  pdbId?: string;
  title?: string;
  classification?: string;
  depositDate?: string;
  method?: string;
  resolution?: number;
  rFactor?: number;
  organism?: string;
  expressionSystem?: string;
  entities?: PdbEntity[];
  citation?: PdbCitation;
  ligands?: PdbLigand[];
  disulfideBondCount?: number;
  atomCount?: number;
  chains?: string[];
}

/** Collect continuation lines for a given record type (e.g. TITLE, EXPDTA).
 *  PDB continuation lines have the same record name with a serial number in cols 9-10. */
function collectLines(lines: string[], recordName: string): string {
  const prefix = recordName.padEnd(6);
  const parts: string[] = [];
  for (const line of lines) {
    if (line.startsWith(prefix))
      parts.push(line.substring(10).trim());
  }
  return parts.join(' ').replace(/\s+/g, ' ').trim();
}

/** Extract a value from SOURCE/COMPND style key-value lines.
 *  These records use "KEY: value;" format across continuation lines. */
function extractKeyValue(lines: string[], recordName: string, key: string): string | undefined {
  const fullText = collectLines(lines, recordName);
  // Match KEY: value; or KEY: value (at end of string)
  const regex = new RegExp(`${key}:\\s*([^;]+)`, 'i');
  const match = fullText.match(regex);
  return match ? match[1].trim() : undefined;
}

/** Parse COMPND records into entities grouped by MOL_ID. */
function parseEntities(lines: string[]): PdbEntity[] | undefined {
  const fullText = collectLines(lines, 'COMPND');
  if (!fullText) return undefined;

  // Split by MOL_ID
  const molBlocks = fullText.split(/MOL_ID:\s*/i).filter(Boolean);
  const entities: PdbEntity[] = [];

  for (const block of molBlocks) {
    const molIdMatch = block.match(/^(\d+)\s*;/);
    const molId = molIdMatch ? molIdMatch[1] : (entities.length + 1).toString();

    const moleculeMatch = block.match(/MOLECULE:\s*([^;]+)/i);
    const chainMatch = block.match(/CHAIN:\s*([^;]+)/i);

    entities.push({
      molId,
      molecule: moleculeMatch ? moleculeMatch[1].trim() : undefined,
      chains: chainMatch ? chainMatch[1].split(',').map((c) => c.trim()).filter(Boolean) : undefined,
    });
  }

  return entities.length > 0 ? entities : undefined;
}

/** Parse JRNL records into citation info. */
function parseCitation(lines: string[]): PdbCitation | undefined {
  const jrnlLines = lines.filter((l) => l.startsWith('JRNL  '));
  if (jrnlLines.length === 0) return undefined;

  const titleParts: string[] = [];
  const refParts: string[] = [];
  let doi: string | undefined;
  let pmid: number | undefined;

  for (const line of jrnlLines) {
    const subtype = line.substring(12, 16).trim();
    const content = line.substring(19).trim();
    if (subtype === 'TITL')
      titleParts.push(content);
    else if (subtype === 'REF')
      refParts.push(content);
    else if (subtype === 'DOI')
      doi = content;
    else if (subtype === 'PMID')
      pmid = parseInt(content, 10) || undefined;
  }

  const title = titleParts.join(' ').replace(/\s+/g, ' ').trim() || undefined;
  const journal = refParts.length > 0 ?
    refParts[0].replace(/\s{2,}/g, ' ').trim().split(/\s{2,}|V\.\s/)[0].trim() :
    undefined;

  if (!title && !journal && !doi && !pmid) return undefined;
  return {title, journal, doi, pmid};
}

/** Parse HET/HETNAM/FORMUL records into ligand info. */
function parseLigands(lines: string[]): PdbLigand[] | undefined {
  const ligandMap: {[id: string]: PdbLigand} = {};

  // HET records: residue name, chain, number of HETATM atoms
  for (const line of lines) {
    if (line.startsWith('HET   ')) {
      const id = line.substring(7, 10).trim();
      const chain = line.substring(12, 13).trim() || undefined;
      if (id && id !== 'HOH') {
        if (!ligandMap[id])
          ligandMap[id] = {id, chain};
      }
    }
  }

  // HETNAM records: compound name
  for (const line of lines) {
    if (line.startsWith('HETNAM')) {
      const id = line.substring(11, 14).trim();
      const namePart = line.substring(15).trim();
      if (id && ligandMap[id]) {
        ligandMap[id].name = ligandMap[id].name ?
          ligandMap[id].name + ' ' + namePart :
          namePart;
      }
    }
  }

  // FORMUL records: chemical formula
  for (const line of lines) {
    if (line.startsWith('FORMUL')) {
      const id = line.substring(12, 15).trim();
      const formulaPart = line.substring(18).trim();
      if (id && ligandMap[id]) {
        ligandMap[id].formula = ligandMap[id].formula ?
          ligandMap[id].formula + ' ' + formulaPart :
          formulaPart;
      }
    }
  }

  const ligands = Object.values(ligandMap);
  return ligands.length > 0 ? ligands : undefined;
}

/** Parse all available header information from a PDB file string. */
export function parsePdbHeaders(pdbText: string): PdbHeaderInfo {
  const lines = pdbText.split('\n');
  const result: PdbHeaderInfo = {};

  // HEADER: classification (cols 11-50), date (cols 51-59), PDB ID (cols 63-66)
  const headerLine = lines.find((l) => l.startsWith('HEADER'));
  if (headerLine) {
    const classification = headerLine.substring(10, 50).trim();
    if (classification) result.classification = classification;

    const dateStr = headerLine.substring(50, 59).trim();
    if (dateStr && /\d{2}-[A-Z]{3}-\d{2}/.test(dateStr)) result.depositDate = dateStr;

    const pdbId = headerLine.substring(62, 66).trim();
    if (pdbId && /^[A-Za-z0-9]{4}$/.test(pdbId)) result.pdbId = pdbId;
  }

  // TITLE
  const title = collectLines(lines, 'TITLE ');
  if (title) result.title = title;

  // EXPDTA
  const method = collectLines(lines, 'EXPDTA');
  if (method) result.method = method;

  // REMARK 2: resolution
  for (const line of lines) {
    if (line.startsWith('REMARK   2') && line.includes('RESOLUTION')) {
      const match = line.match(/(\d+\.\d+)\s+ANGSTROMS/);
      if (match) result.resolution = parseFloat(match[1]);
      break;
    }
  }

  // REMARK 3: R-factor
  for (const line of lines) {
    if (line.startsWith('REMARK   3') && /R VALUE\s+\(WORKING SET\)/i.test(line)) {
      const match = line.match(/:\s*(\d+\.\d+)/);
      if (match) result.rFactor = parseFloat(match[1]);
      break;
    }
  }

  // SOURCE: organism, expression system
  const organism = extractKeyValue(lines, 'SOURCE', 'ORGANISM_SCIENTIFIC');
  if (organism) result.organism = organism;

  const exprSystem = extractKeyValue(lines, 'SOURCE', 'EXPRESSION_SYSTEM(?!_)');
  // More precise: get the value directly
  const sourceText = collectLines(lines, 'SOURCE');
  if (sourceText) {
    const exprMatch = sourceText.match(/EXPRESSION_SYSTEM:\s*([^;]+)/i);
    if (exprMatch) {
      // Make sure we don't capture EXPRESSION_SYSTEM_TAXID etc.
      result.expressionSystem = exprMatch[1].trim();
    }
  }

  // COMPND: entities
  result.entities = parseEntities(lines);

  // JRNL: citation
  result.citation = parseCitation(lines);

  // Ligands: HET, HETNAM, FORMUL
  result.ligands = parseLigands(lines);

  // SSBOND count
  const ssbondCount = lines.filter((l) => l.startsWith('SSBOND')).length;
  if (ssbondCount > 0) result.disulfideBondCount = ssbondCount;

  // Atom count and chains from ATOM/HETATM records
  let atomCount = 0;
  const chainSet = new Set<string>();
  for (const line of lines) {
    if (line.startsWith('ATOM  ') || line.startsWith('HETATM')) {
      atomCount++;
      const chain = line.substring(21, 22).trim();
      if (chain) chainSet.add(chain);
    }
  }
  if (atomCount > 0) result.atomCount = atomCount;
  if (chainSet.size > 0) result.chains = [...chainSet].sort();

  return result;
}
