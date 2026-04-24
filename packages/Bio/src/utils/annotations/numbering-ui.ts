/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS, ALIGNMENT, ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  SeqAnnotation, SeqAnnotationHit, AnnotationCategory,
} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
import {
  setColumnAnnotations, getColumnAnnotations,
  getOrCreateAnnotationColumn, getRowAnnotations, setRowAnnotations, mergeRowHits,
} from './annotation-manager';
/** An engine entry — a dynamically discovered DG.Func with meta.role: 'antibodyNumbering'. */
interface NumberingEngine {
  label: string;
  /** Unique key — `${package}:${name}`. */
  key: string;
  func: DG.Func;
}

/** Discovers all registered antibody numbering engines (functions with
 *  meta.role = 'antibodyNumbering'). Built-in engines: Bio package ships an
 *  immunum-WASM-based engine; other packages (Proteomics etc.) can register
 *  AntPack/ANARCI/etc. by adding their own meta.role function. */
function discoverEngines(): NumberingEngine[] {
  const engines: NumberingEngine[] = [];
  const funcs = DG.Func.find({meta: {role: 'antibodyNumbering'}});
  if (funcs.length === 0) {
    grok.shell.error('No antibody numbering engines found. Make sure the Bio package is up to date.');
    throw new Error('No antibody numbering engines found.');
  }
  for (const f of funcs) {
    const pkgName = f.package?.name ?? '';
    engines.push({
      label: f.friendlyName || f.name,
      key: pkgName ? `${pkgName}:${f.name}` : f.name,
      func: f,
    });
  }
  return engines;
}

/** Reads the `scheme` parameter's declared choices off an engine function.
 *  Engines advertise which schemes they support via the `choices` option of
 *  their `scheme: string` parameter; if the function has none we fall back to
 *  a single generic IMGT option so the dropdown is never empty. */
function getEngineSchemes(engine: NumberingEngine): string[] {
  const schemeInput = engine.func.inputs.find((p) => p.name.toLowerCase() === 'scheme');
  const choices = schemeInput?.choices ?? [];
  return choices.length > 0 ? choices.slice() : ['imgt'];
}

export function showNumberingSchemeDialog(): void {
  const df = grok.shell.tv?.dataFrame;
  if (!df) {
    grok.shell.warning('No table open');
    return;
  }

  const seqCols = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (seqCols.length === 0) {
    grok.shell.warning('No macromolecule columns found');
    return;
  }

  const engines = discoverEngines();
  const engineLabels = engines.map((e) => e.label);

  const tableInput = ui.input.table('Table', {value: df});
  const seqInput = ui.input.column('Sequence', {
    table: df, value: seqCols[0],
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
  });
  const engineInput = ui.input.choice('Engine', {
    value: engineLabels[0], items: engineLabels,
  });
  const initialSchemes = getEngineSchemes(engines[0]);
  const schemeInput = ui.input.choice('Scheme', {value: initialSchemes[0], items: initialSchemes});

  // Switch the scheme list when the user picks a different engine — each engine
  // advertises its own set via the `scheme` parameter's `choices`.
  engineInput.onChanged.subscribe(() => {
    const selected = engines.find((e) => e.label === engineInput.value);
    if (!selected) return;
    const schemes = getEngineSchemes(selected);
    const prev = schemeInput.value;
    schemeInput.items = schemes;
    schemeInput.value = schemes.includes(prev ?? '') ? prev : schemes[0];
  });

  const dialog = ui.dialog({title: 'Apply Antibody Numbering'})
    .add(ui.inputs([tableInput, seqInput, engineInput, schemeInput]))
    .onOK(async () => {
      const seqCol = seqInput.value!;
      const schemeName = schemeInput.value!;
      const selectedLabel = engineInput.value!;
      const engine = engines.find((e) => e.label === selectedLabel) ?? engines[0];
      const pi = DG.TaskBarProgressIndicator.create(`Applying ${schemeName} numbering...`);
      try {
        const result: DG.DataFrame = await engine.func.apply({df: df, seqCol: seqCol, scheme: schemeName});
        applyNumberingResults(df, seqCol, result, schemeName, true, engine.label);
      } catch (err: any) {
        grok.shell.error(`Numbering failed: ${err.message ?? err}`);
        console.error(err);
      } finally {
        pi.close();
      }
    });

  dialog.show();
}

/** Builds a map from ungapped character index to gapped character index.
 *  Used when the source column is already aligned (MSA) — numbering engines strip gaps,
 *  so their output indices refer to the ungapped sequence, not the gapped original. */
function buildUngappedToGappedMap(gappedSeq: string): number[] {
  const map: number[] = [];
  for (let g = 0; g < gappedSeq.length; g++) {
    if (gappedSeq[g] !== '-' && gappedSeq[g] !== '.')
      map.push(g);
  }
  return map;
}

/** Parses a position code into [numericPart, insertionLetter] for sorting.
 *  E.g. "27" → [27, ""], "111A" → [111, "A"], "27B" → [27, "B"]. */
function parsePositionCode(code: string): [number, string] {
  const match = code.match(/^(\d+)([A-Z]?)$/);
  if (!match) return [Infinity, code];
  return [parseInt(match[1], 10), match[2]];
}

/** Sorts position codes in scheme order: numeric ascending, then insertion letter. */
function sortPositionCodes(codes: string[]): string[] {
  return codes.slice().sort((a, b) => {
    const [numA, insA] = parsePositionCode(a);
    const [numB, insB] = parsePositionCode(b);
    if (numA !== numB) return numA - numB;
    return insA.localeCompare(insB);
  });
}

/** Builds unified position list from all rows and creates an aligned sequence column.
 *  Includes flanking residues (before/after the numbered region) padded with gaps.
 *  Layout: [pre-region gaps+residues] [scheme-aligned region] [post-region residues+gaps]
 *  @returns aligned column, full position list (including flanking), and the pre-region offset. */
function createAlignedColumn(
  df: DG.DataFrame, seqCol: DG.Column<string>, result: DG.DataFrame,
): {alignedCol: DG.Column<string>; unifiedPositions: string[]; preOffset: number} | null {
  const numberingMapCol = result.col('numbering_map');
  if (!numberingMapCol) return null;

  // Pass 1: collect all scheme position codes and per-row flanking lengths
  const allCodes = new Set<string>();
  const rowMaps: (Record<string, number> | null)[] = [];
  const rowPreLens: number[] = []; // chars before first numbered position
  const rowPostLens: number[] = []; // chars after last numbered position

  for (let i = 0; i < result.rowCount; i++) {
    const mapStr = numberingMapCol.get(i);
    const rawSeq = seqCol.get(i) ?? '';
    if (!mapStr) {
      rowMaps.push(null);
      rowPreLens.push(0);
      rowPostLens.push(0);
      continue;
    }
    try {
      const posToCharIdx: Record<string, number> = JSON.parse(mapStr);
      rowMaps.push(posToCharIdx);
      for (const code of Object.keys(posToCharIdx))
        allCodes.add(code);

      // Find min/max char indices that were numbered
      const charIndices = Object.values(posToCharIdx);
      const minChar = Math.min(...charIndices);
      const maxChar = Math.max(...charIndices);
      rowPreLens.push(minChar); // chars before first numbered
      rowPostLens.push(Math.max(0, rawSeq.length - maxChar - 1)); // chars after last numbered
    } catch {
      rowMaps.push(null);
      rowPreLens.push(0);
      rowPostLens.push(0);
    }
  }

  if (allCodes.size === 0) return null;

  const maxPreLen = Math.max(0, ...rowPreLens);
  const maxPostLen = Math.max(0, ...rowPostLens);

  // Build position names: [pre-flanking] + [scheme positions] + [post-flanking]
  const schemePositions = sortPositionCodes(Array.from(allCodes));
  const preNames: string[] = [];
  for (let p = maxPreLen; p > 0; p--)
    preNames.push(`N-${p}`);
  const postNames: string[] = [];
  for (let p = 1; p <= maxPostLen; p++)
    postNames.push(`C+${p}`);

  const unifiedPositions = [...preNames, ...schemePositions, ...postNames];
  const totalLen = unifiedPositions.length;
  const preOffset = maxPreLen; // scheme region starts at this index in the aligned string

  // Map scheme position codes → index in the full unified list
  const posToUnifiedIdx = new Map<string, number>();
  for (let s = 0; s < schemePositions.length; s++)
    posToUnifiedIdx.set(schemePositions[s], preOffset + s);

  // Pass 2: build aligned sequences
  const colName = df.columns.getUnusedName(`${seqCol.name} (aligned)`);
  const alignedCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, colName, df.rowCount);

  for (let i = 0; i < df.rowCount; i++) {
    const map = i < rowMaps.length ? rowMaps[i] : null;
    const rawSeq = seqCol.get(i) ?? '';

    if (!map) {
      alignedCol.set(i, '-'.repeat(totalLen));
      continue;
    }

    const aligned = new Array<string>(totalLen).fill('-');

    // Place scheme-numbered residues
    for (const [posCode, charIdx] of Object.entries(map)) {
      const uIdx = posToUnifiedIdx.get(posCode);
      if (uIdx != null && charIdx < rawSeq.length)
        aligned[uIdx] = rawSeq[charIdx];
    }

    // Place pre-region flanking residues (right-aligned within the pre-region block)
    const preLen = rowPreLens[i];
    const charIndices = Object.values(map);
    const minChar = Math.min(...charIndices);
    for (let p = 0; p < preLen; p++)
      aligned[preOffset - preLen + p] = rawSeq[minChar - preLen + p];

    // Place post-region flanking residues (left-aligned within the post-region block)
    const postLen = rowPostLens[i];
    const maxChar = Math.max(...charIndices);
    const postStart = preOffset + schemePositions.length;
    for (let p = 0; p < postLen; p++)
      aligned[postStart + p] = rawSeq[maxChar + 1 + p];

    alignedCol.set(i, aligned.join(''));
  }

  // Set macromolecule tags on the aligned column
  alignedCol.semType = DG.SEMTYPE.MACROMOLECULE;
  alignedCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  alignedCol.setTag(bioTAGS.alphabet, ALPHABET.PT);
  alignedCol.meta.units = NOTATION.FASTA;
  alignedCol.setTag(DG.Tags.CellRenderer, 'sequence');
  alignedCol.setTag(bioTAGS.positionNames, unifiedPositions.join(', '));

  return {alignedCol, unifiedPositions, preOffset};
}

/** Applies numbering results (from either engine) to the sequence column and dataframe.
 *
 *  Annotation strategy:
 *  - Original column: row-level region spans only (char indices from numbering_map).
 *    No column-level position names (they differ per row).
 *  - Aligned column: column-level annotations only (all rows share unified positions).
 *    Position names tag set to the unified list. */
function applyNumberingResults(
  df: DG.DataFrame, seqCol: DG.Column<string>, result: DG.DataFrame,
  schemeName: string, populateRegions: boolean, engineLabel: string,
): void {
  if (!result || result.rowCount === 0) {
    grok.shell.warning('No numbering results returned');
    return;
  }

  const posNamesCol = result.getCol('position_names');
  const chainTypeCol = result.getCol('chain_type');
  const annotJsonCol = result.getCol('annotations_json');
  const numberingMapCol = result.col('numbering_map');

  // If the source column is already aligned (MSA), remap numbering_map indices
  // from ungapped to gapped, since numbering engines strip gaps before processing.
  const isAligned = true; // always treat as aligned to handle remapping
  if (isAligned && numberingMapCol) {
    for (let i = 0; i < result.rowCount; i++) {
      const mapStr = numberingMapCol.get(i);
      if (!mapStr) continue;
      try {
        const posToCharIdx: Record<string, number> = JSON.parse(mapStr);
        const rawSeq = seqCol.get(i) ?? '';
        const ungapToGap = buildUngappedToGappedMap(rawSeq);
        const remapped: Record<string, number> = {};
        for (const [posCode, ungappedIdx] of Object.entries(posToCharIdx)) {
          if (ungappedIdx < ungapToGap.length)
            remapped[posCode] = ungapToGap[ungappedIdx];
        }
        numberingMapCol.set(i, JSON.stringify(remapped));
      } catch { /* skip */ }
    }
  }

  // Pick the row with the most annotations for column-level representative data
  let bestRowIdx = -1;
  let bestAnnotCount = -1;

  for (let i = 0; i < result.rowCount; i++) {
    const pn = posNamesCol.get(i);
    if (!pn || pn.length === 0) continue;
    const aj = annotJsonCol.get(i);
    let count = 0;
    if (aj)
      try { count = JSON.parse(aj).length; } catch { /* skip */ }
    if (count > bestAnnotCount) {
      bestAnnotCount = count;
      bestRowIdx = i;
    }
  }

  const chainType = bestRowIdx >= 0 ? (chainTypeCol.get(bestRowIdx) ?? '') : '';
  const annotationsJson = bestRowIdx >= 0 ? (annotJsonCol.get(bestRowIdx) ?? '[]') : '[]';

  if (bestRowIdx < 0) {
    grok.shell.warning(`${engineLabel} could not number the sequences. Check that they are valid antibody variable region sequences.`);
    return;
  }

  // Mark scheme on original column (no position names — they differ per row)
  seqCol.setTag(bioTAGS.numberingScheme, schemeName);

  // --- Original column: column-level annotation definitions (needed for renderer
  // to resolve annotation IDs → colors/names) + row-level region spans ---
  if (populateRegions) {
    try {
      const regionAnnotations: SeqAnnotation[] = JSON.parse(annotationsJson);
      const existing = getColumnAnnotations(seqCol).filter((a) => a.category !== AnnotationCategory.Structure);
      setColumnAnnotations(seqCol, [...existing, ...regionAnnotations]);
    } catch (err) {
      console.warn('Failed to set annotation definitions on original column:', err);
    }
  }

  if (populateRegions && numberingMapCol) {
    try {
      const annotCol = getOrCreateAnnotationColumn(df, seqCol);
      for (let i = 0; i < result.rowCount; i++) {
        const mapStr = numberingMapCol.get(i);
        const rowAnnotJson = annotJsonCol.get(i);
        if (!mapStr || !rowAnnotJson) continue;
        const posToCharIdx: Record<string, number> = JSON.parse(mapStr);
        const rowRegions: SeqAnnotation[] = JSON.parse(rowAnnotJson);

        const regionHits: SeqAnnotationHit[] = [];
        for (const region of rowRegions) {
          if (region.start == null || region.end == null) continue;
          const startCharIdx = posToCharIdx[region.start];
          const endCharIdx = posToCharIdx[region.end];
          if (startCharIdx == null || endCharIdx == null) continue;
          regionHits.push({
            annotationId: region.id,
            positionIndex: startCharIdx,
            endPositionIndex: endCharIdx,
            positionName: region.start,
            matchedMonomers: '',
          });
        }

        const existingHits = getRowAnnotations(annotCol, i) ?? [];
        setRowAnnotations(annotCol, i, mergeRowHits(existingHits, regionHits, true, false));
      }
    } catch (err) {
      console.warn('Failed to store per-row region data on original column:', err);
    }
  }

  // --- Aligned column: column-level annotations only ---
  const alignment = createAlignedColumn(df, seqCol, result);
  if (alignment) {
    df.columns.insert(alignment.alignedCol, df.columns.toList().indexOf(seqCol) + 1);
    if (grok.shell.tv?.dataFrame === df)
      grok.shell.tv.grid.scrollToCell(seqCol, 0);

    alignment.alignedCol.setTag(bioTAGS.numberingScheme, schemeName);

    if (populateRegions) {
      try {
        const regionAnnotations: SeqAnnotation[] = JSON.parse(annotationsJson);
        const existing = getColumnAnnotations(alignment.alignedCol).filter((a) => a.category !== AnnotationCategory.Structure);
        setColumnAnnotations(alignment.alignedCol, [...existing, ...regionAnnotations]);
        // chunk for vd regions viewer if that becomes a desired feature in the future
        // if (grok.shell.tv?.dataFrame === df) {
        //   (async () => {
        //     const vdRegionsViewer: VdRegionsViewer = await grok.shell.tv.dataFrame.plot.fromType('VdRegions',) as VdRegionsViewer;
        //     vdRegionsViewer.chains = [chainType];
        //     vdRegionsViewer.setData(regionAnnotations.map((a, i) => ({
        //       type: a.name?.toLowerCase().startsWith('cdr') ? VdRegionType.CDR : VdRegionType.FR,
        //       name: a.name,
        //       chain: chainType,
        //       order: i,
        //       sequenceColumnName: alignment.alignedCol.name,
        //       positionStartName: a.start ?? '',
        //       positionEndName: a.end ?? '',
        //     } satisfies VdRegion)));
        //     grok.shell.tv.addViewer(vdRegionsViewer);
        //   })();
        // }
      } catch (err) {
        console.warn('Failed to set column-level annotations on aligned column:', err);
      }
    }
  }

  df.fireValuesChanged();

  grok.shell.info(`Numbering applied: ${schemeName}, chain type: ${chainType}`);
}
