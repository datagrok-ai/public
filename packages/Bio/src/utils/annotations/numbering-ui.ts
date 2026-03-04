/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  SeqAnnotation, SeqAnnotationHit, AnnotationCategory,
} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
import {NumberingScheme} from '@datagrok-libraries/bio/src/utils/macromolecule/numbering-schemes';
import {
  setColumnAnnotations, getColumnAnnotations,
  getOrCreateAnnotationColumn, getRowAnnotations, setRowAnnotations, mergeRowHits,
} from './annotation-manager';
import {_package} from '../../package';
import type {NumberingResult, Scheme} from '../antibody-numbering (WIP)';

const BUILTIN_ENGINE_KEY = '__builtin__';
const BUILTIN_ENGINE_LABEL = 'Built-in (TypeScript)';

/** An engine entry: either a dynamically discovered DG.Func or the built-in TS engine. */
interface NumberingEngine {
  /** Display label for the dropdown */
  label: string;
  /** Unique key — nqName for DG.Func engines, BUILTIN_ENGINE_KEY for built-in */
  key: string;
  /** The DG.Func to call, or null for the built-in engine */
  func: DG.Func | null;
}

/** Discovers all registered antibody numbering engines + the built-in TS engine.
 *  Dynamic engines (meta.role = 'antibodyNumbering') come first; built-in is last. */
function discoverEngines(): NumberingEngine[] {
  const engines: NumberingEngine[] = [];

  const funcs = DG.Func.find({meta: {role: 'antibodyNumbering'}});
  for (const f of funcs) {
    const pkgName = f.package?.name ?? '';
    const label = f.friendlyName || f.name;
    engines.push({
      label: pkgName ? `${label} (${pkgName})` : label,
      key: pkgName ? `${pkgName}:${f.name}` : f.name,
      func: f,
    });
  }

  // Built-in TS engine is always last
  engines.push({label: BUILTIN_ENGINE_LABEL, key: BUILTIN_ENGINE_KEY, func: null});
  return engines;
}

/** Converts TS NumberingResult[] to a DG.DataFrame matching the expected output shape.
 *  Columns: position_names, chain_type, annotations_json, numbering_detail, numbering_map. */
export function numberingResultsToDataFrame(results: NumberingResult[]): DG.DataFrame {
  const n = results.length;
  const posNames = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'position_names', n);
  const chainTypes = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'chain_type', n);
  const annotJson = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'annotations_json', n);
  const numDetail = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'numbering_detail', n);
  const numMap = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'numbering_map', n);

  for (let i = 0; i < n; i++) {
    const r = results[i];
    if (r.error && r.percentIdentity < 0.3) {
      posNames.set(i, '');
      chainTypes.set(i, '');
      annotJson.set(i, '[]');
      numDetail.set(i, '');
      numMap.set(i, '');
    } else {
      posNames.set(i, r.positionNames);
      chainTypes.set(i, r.chainType);
      annotJson.set(i, JSON.stringify(r.annotations));
      numDetail.set(i, JSON.stringify(r.numberingDetail));
      numMap.set(i, JSON.stringify(r.numberingMap));
    }
  }

  return DG.DataFrame.fromColumns([posNames, chainTypes, annotJson, numDetail, numMap]);
}

/** Runs the built-in TS numbering engine on all rows of a sequence column. */
async function runBuiltinNumbering(
  seqCol: DG.Column<string>, schemeName: string,
): Promise<DG.DataFrame> {
  const {numberSequences, extractSequence} = await import('../antibody-numbering (WIP)');
  const scheme = schemeName.toLowerCase() as Scheme;

  const sequences: string[] = [];
  for (let i = 0; i < seqCol.length; i++) {
    const raw = seqCol.get(i);
    sequences.push(extractSequence(raw ?? ''));
  }

  const results = numberSequences(sequences, scheme);
  return numberingResultsToDataFrame(results);
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
  const schemeChoices = Object.values(NumberingScheme);

  const tableInput = ui.input.table('Table', {value: df});
  const seqInput = ui.input.column('Sequence', {
    table: df, value: seqCols[0],
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
  });
  const schemeInput = ui.input.choice('Scheme', {value: NumberingScheme.IMGT, items: schemeChoices});
  const engineInput = ui.input.choice('Engine', {
    value: engineLabels[0], items: engineLabels,
  });
  // const populateRegions = ui.input.bool('Populate FR/CDR regions', {value: true});
  // const openVdRegions = ui.input.bool('Open VD Regions viewer', {value: true});

  const dialog = ui.dialog({title: 'Apply Antibody Numbering'})
    .add(ui.inputs([tableInput, seqInput, schemeInput, engineInput]))
    .onOK(async () => {
      const seqCol = seqInput.value!;
      const schemeName = schemeInput.value!;
      const selectedLabel = engineInput.value!;
      const engine = engines.find((e) => e.label === selectedLabel) ?? engines[engines.length - 1];
      const pi = DG.TaskBarProgressIndicator.create(`Applying ${schemeName} numbering...`);
      try {
        let result: DG.DataFrame;
        if (engine.func)
          result = await engine.func.apply({df: df, seqCol: seqCol, scheme: schemeName.toLowerCase()});
        else
          result = await runBuiltinNumbering(seqCol, schemeName);

        applyNumberingResults(df, seqCol, result, schemeName, true, engine.label);

        // // Open VD Regions viewer
        // if (openVdRegions.value && grok.shell.tv) {
        //   try {
        //     await grok.shell.tv.dataFrame.plot.fromType('VdRegions', {});
        //   } catch (err) {
        //     console.warn('Could not open VD Regions viewer:', err);
        //   }
        // }
      } catch (err: any) {
        grok.shell.error(`Numbering failed: ${err.message ?? err}`);
        console.error(err);
      } finally {
        pi.close();
      }
    });

  dialog.show();
}

/** Applies numbering results (from either engine) to the sequence column and dataframe. */
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

  // Pick the row with the most annotations for column-level data.
  // Some rows may only have partial numbering (e.g. FR1+CDR1 only),
  // so we want the most complete one as the column-level representative.
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
      break; // Prefer the first complete one, no need to scan all rows for now. remove if necessary
    }
  }

  const positionNames = bestRowIdx >= 0 ? (posNamesCol.get(bestRowIdx) ?? '') : '';
  const chainType = bestRowIdx >= 0 ? (chainTypeCol.get(bestRowIdx) ?? '') : '';
  const annotationsJson = bestRowIdx >= 0 ? (annotJsonCol.get(bestRowIdx) ?? '[]') : '[]';

  if (!positionNames) {
    grok.shell.warning(`${engineLabel} could not number the sequences. Check that they are valid antibody variable region sequences.`);
    return;
  }

  // Set position names
  seqCol.setTag(bioTAGS.positionNames, positionNames);
  seqCol.setTag(bioTAGS.numberingScheme, schemeName);

  // Set column-level region annotations
  if (populateRegions) {
    try {
      const regionAnnotations: SeqAnnotation[] = JSON.parse(annotationsJson);
      const existing = getColumnAnnotations(seqCol).filter((a) => a.category !== AnnotationCategory.Structure);
      setColumnAnnotations(seqCol, [...existing, ...regionAnnotations]);
    } catch (err) {
      console.warn('Failed to parse region annotations:', err);
    }
  }

  // Store per-row region spans in the companion annotation column
  if (populateRegions && numberingMapCol) {
    try {
      const annotCol = getOrCreateAnnotationColumn(df, seqCol);
      for (let i = 0; i < result.rowCount; i++) {
        const mapStr = numberingMapCol.get(i);
        const rowAnnotJson = annotJsonCol.get(i);
        if (!mapStr || !rowAnnotJson) continue;
        const posToCharIdx: Record<string, number> = JSON.parse(mapStr);
        const rowRegions: SeqAnnotation[] = JSON.parse(rowAnnotJson);

        // Convert region annotations to row-level hits with character indices
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

        // Merge: preserve existing liability hits, replace region hits
        const existingHits = getRowAnnotations(annotCol, i) ?? [];
        setRowAnnotations(annotCol, i, mergeRowHits(existingHits, regionHits, true, false));
      }
    } catch (err) {
      console.warn('Failed to store per-row region data:', err);
    }
  }

  df.fireValuesChanged();
  grok.shell.info(`Numbering applied: ${schemeName}, chain type: ${chainType}`);
}
