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

const ENGINE_BUILTIN = 'Built-in (TypeScript)';
const ENGINE_ANTPACK = 'AntPack (Python)';

/** Converts TS NumberingResult[] to a DG.DataFrame matching the Python script output shape.
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

  const schemeChoices = Object.values(NumberingScheme);

  const tableInput = ui.input.table('Table', {value: df});
  const seqInput = ui.input.column('Sequence', {
    table: df, value: seqCols[0],
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
  });
  const schemeInput = ui.input.choice('Scheme', {value: NumberingScheme.IMGT, items: schemeChoices});
  const engineInput = ui.input.choice('Engine', {
    value: ENGINE_BUILTIN, items: [ENGINE_BUILTIN, ENGINE_ANTPACK],
  });
  const populateRegions = ui.input.bool('Populate FR/CDR regions', {value: true});
  const openVdRegions = ui.input.bool('Open VD Regions viewer', {value: true});

  const hintDiv = ui.div([
    ui.divText('Built-in engine runs in-browser. AntPack requires Python environment.', {style: {fontSize: '11px', color: '#888', marginTop: '8px'}}),
  ]);

  const dialog = ui.dialog({title: 'Apply Antibody Numbering'})
    .add(ui.inputs([tableInput, seqInput, schemeInput, engineInput, populateRegions, openVdRegions]))
    .add(hintDiv)
    .onOK(async () => {
      const seqCol = seqInput.value!;
      const schemeName = schemeInput.value!;
      const engine = engineInput.value!;
      const pi = DG.TaskBarProgressIndicator.create(`Applying ${schemeName} numbering...`);
      try {
        let result: DG.DataFrame;
        if (engine === ENGINE_ANTPACK) {
          result = await grok.functions.call(
            'Bio:NumberAntibodySequences', {
              df: df,
              seqCol: seqCol,
              scheme: schemeName.toLowerCase(),
            },
          );
        } else {
          result = await runBuiltinNumbering(seqCol, schemeName);
        }

        applyNumberingResults(df, seqCol, result, schemeName, populateRegions.value!, engine);

        // Open VD Regions viewer
        if (openVdRegions.value && grok.shell.tv) {
          try {
            await grok.shell.tv.dataFrame.plot.fromType('VdRegions', {});
          } catch (err) {
            console.warn('Could not open VD Regions viewer:', err);
          }
        }
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
  schemeName: string, populateRegions: boolean, engine: string,
): void {
  if (!result || result.rowCount === 0) {
    grok.shell.warning('No numbering results returned');
    return;
  }

  const posNamesCol = result.getCol('position_names');
  const chainTypeCol = result.getCol('chain_type');
  const annotJsonCol = result.getCol('annotations_json');
  const numberingMapCol = result.col('numbering_map');

  // Use the first non-empty result for column-level position names
  let positionNames = '';
  let chainType = '';
  let annotationsJson = '[]';

  for (let i = 0; i < result.rowCount; i++) {
    const pn = posNamesCol.get(i);
    if (pn && pn.length > 0) {
      positionNames = pn;
      chainType = chainTypeCol.get(i) ?? '';
      annotationsJson = annotJsonCol.get(i) ?? '[]';
      break;
    }
  }

  if (!positionNames) {
    const engineLabel = engine === ENGINE_ANTPACK ? 'AntPack' : 'Built-in engine';
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
