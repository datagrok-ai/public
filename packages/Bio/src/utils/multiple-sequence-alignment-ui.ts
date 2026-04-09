/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Subscription} from 'rxjs';

import {ColumnInputOptions} from '@datagrok-libraries/utils/src/type-declarations';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {NotationProviderBase} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {MsaWarning, runKalign, checkForSingleSeqClusters} from './multiple-sequence-alignment';
import {checkInputColumn} from './check-input-column';
import {MultipleSequenceAlignmentUIOptions} from './types';
import {kalignVersion, MSA_ENGINE_ROLE} from './constants';
import {_package} from '../package';

import '../../css/msa.css';

type AlignmentMode = 'kalign' | 'engine';

/** State holder for the MSA dialog, avoids TypeScript narrowing issues with closures. */
class MsaDialogState {
  mode: AlignmentMode = 'kalign';
  currentFunc: DG.Func | null = null;
  currentFuncCall: DG.FuncCall | null = null;
}


export async function multipleSequenceAlignmentUI(
  options: MultipleSequenceAlignmentUIOptions, seqHelper: ISeqHelper,
): Promise<DG.Column> {
  return new Promise(async (resolve, reject) => {
    try {
      const table = options.col?.dataFrame ?? grok.shell.t;
      if (!table) {
        reject(new MsaWarning(ui.divText('MSA requires a dataset with a macromolecule column.')));
        return;
      }

      const seqCol = options.col ?? table.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
      if (!seqCol) {
        reject(new MsaWarning(ui.divText('MSA requires a dataset with a macromolecule column.')));
        return;
      }

      const state = new MsaDialogState();

      // --- Common UI ---

      let prevSeqCol = seqCol;
      const colInput = ui.input.column('Sequence', {
        table, value: seqCol,
        onValueChanged: async (value: DG.Column<any>) => {
          if (!value || value.semType !== DG.SEMTYPE.MACROMOLECULE) {
            okBtn.disabled = true;
            await DG.delay(0);
            colInput.value = prevSeqCol as DG.Column<string>;
            return;
          }
          prevSeqCol = value;
          okBtn.disabled = false;
          await onColumnChanged(value);
        },
        filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
      } as ColumnInputOptions) as DG.InputBase<DG.Column<string>>;
      colInput.setTooltip('Sequences column to use for alignment');

      const clustersColInput = ui.input.column('Clusters', {table, value: options.clustersCol!});
      clustersColInput.nullable = true;

      const onlySelectedInput = ui.input.bool('Selected Rows Only', {value: false});

      // --- Kalign UI (canonical sequences) ---

      const kalignGapOpen = ui.input.float('Gap open', {value: options?.kalign?.gapOpen});
      kalignGapOpen.setTooltip('Gap opening penalty at group-to-group alignment');
      const kalignGapExtend = ui.input.float('Gap extend', {value: options?.kalign?.gapExtend});
      kalignGapExtend.setTooltip('Gap extension penalty to skip the alignment');
      const kalignTerminalGap = ui.input.float('Terminal gap', {value: options?.kalign?.terminalGap});
      kalignTerminalGap.setTooltip('Penalty for opening a gap at the beginning or end of the sequence');
      const kalignVersionDiv = ui.p(`Kalign version: ${kalignVersion}`, 'kalign-version');

      const kalignParamsDiv = ui.inputs([kalignGapOpen, kalignGapExtend, kalignTerminalGap]);
      kalignParamsDiv.hidden = true;
      const kalignParamsButton = ui.button('Alignment parameters', () => {
        kalignParamsDiv.hidden = !kalignParamsDiv.hidden;
        [kalignGapOpen, kalignGapExtend, kalignTerminalGap].forEach((input) => {
          input.root.style.removeProperty('max-width');
          input.captionLabel.style.removeProperty('max-width');
        });
      }, 'Adjust alignment parameters such as penalties for opening and extending gaps');
      kalignParamsButton.classList.add('msa-params-button');
      kalignParamsButton.prepend(ui.icons.settings(() => null));

      const kalignElements = [kalignParamsDiv, kalignParamsButton, kalignVersionDiv];

      // --- Engine UI (non-canonical sequences, dynamically discovered) ---

      const msaEngines = DG.Func.find({meta: {role: MSA_ENGINE_ROLE}});
      // Sort so the default engine (meta.defaultAlignment) comes first
      msaEngines.sort((a, b) => {
        const aDefault = a.options['defaultAlignment'] === 'true' ? 1 : 0;
        const bDefault = b.options['defaultAlignment'] === 'true' ? 1 : 0;
        return bDefault - aDefault;
      });

      const engineInput = ui.input.choice('Engine', {
        value: msaEngines.length > 0 ? msaEngines[0].friendlyName : '',
        items: msaEngines.map((f) => f.friendlyName),
      });

      const engineDescDiv = ui.div('', {style: {fontSize: '12px', opacity: '0.7', marginBottom: '6px'}});
      const engineParamsDiv = ui.div();
      const engineParamsButton = ui.button('Alignment parameters', () => {
        engineParamsDiv.hidden = !engineParamsDiv.hidden;
      }, 'Adjust engine-specific alignment parameters');
      engineParamsButton.classList.add('msa-params-button');
      engineParamsButton.prepend(ui.icons.settings(() => null));

      // "Include HELM" checkbox: shown when column has a notation provider with fromHelm
      const includeHelmInput = ui.input.bool('Include HELM', {value: true});
      includeHelmInput.setTooltip('Also add the aligned HELM column alongside the converted notation');
      includeHelmInput.root.style.display = 'none';

      const engineElements = [
        engineInput.root, engineDescDiv, engineParamsButton, engineParamsDiv, includeHelmInput.root,
      ];

      // --- State management ---

      async function updateEngineEditor(): Promise<void> {
        engineParamsDiv.innerHTML = '';
        state.currentFuncCall = null;

        const selectedName = engineInput.value;
        state.currentFunc = msaEngines.find((f) => f.friendlyName === selectedName) ?? null;
        engineDescDiv.textContent = state.currentFunc?.description ?? '';
        if (!state.currentFunc) return;

        state.currentFuncCall = state.currentFunc.prepare({});
        const inputs = await state.currentFuncCall.buildEditor(engineParamsDiv);

        // Hide the first input (sequence column) - managed by the dialog's column selector
        if (inputs.length > 0 && inputs[0].inputType === 'column')
          inputs[0].root.style.display = 'none';
      }

      const _engineSub: Subscription = engineInput.onChanged.subscribe(() => updateEngineEditor());

      function switchMode(newMode: AlignmentMode): void {
        state.mode = newMode;
        for (const el of kalignElements)
          el.style.display = newMode === 'kalign' ? '' : 'none';
        for (const el of engineElements)
          el.style.display = newMode === 'engine' ? '' : 'none';
      }

      async function onColumnChanged(col: DG.Column<string>): Promise<void> {
        try {
          if (col.semType !== DG.SEMTYPE.MACROMOLECULE) return;

          const isCanonical = checkInputColumn(
            col, col.name, seqHelper,
            [NOTATION.FASTA, NOTATION.SEPARATOR], [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT],
          )[0];
          const isHelm = checkInputColumn(col, col.name, seqHelper, [NOTATION.HELM], [])[0];
          const isSepUnknown = checkInputColumn(
            col, col.name, seqHelper, [NOTATION.SEPARATOR, NOTATION.CUSTOM, NOTATION.BILN], [ALPHABET.UN],
          )[0];

          if (isCanonical) {
            switchMode('kalign');
            kalignGapOpen.value = null;
            kalignGapExtend.value = null;
            kalignTerminalGap.value = null;
          } else if (isHelm || isSepUnknown) {
            if (msaEngines.length === 0) {
              grok.shell.warning('No MSA engines found for non-canonical sequences.');
              switchMode('kalign');
              return;
            }
            switchMode('engine');
            await updateEngineEditor();

            // Show "Include HELM" checkbox if the column has a notation provider with fromHelm
            const np = col.temp?.[SeqTemps.notationProvider];
            const npCons = np ? np.constructor as typeof NotationProviderBase : null;
            const hasFromHelm = npCons?.implementsFromHelm === true;
            includeHelmInput.root.style.display = hasFromHelm ? '' : 'none';
          } else
            switchMode('kalign');
        } catch (err: any) {
          const errMsg = err instanceof Error ? err.message : err.toString();
          grok.shell.error(errMsg);
          _package.logger.error(errMsg);
        }
      }

      // --- Alignment execution ---

      async function doAlignment(): Promise<DG.Column<string>> {
        const col = colInput.value;
        if (!col || col.semType !== DG.SEMTYPE.MACROMOLECULE)
          throw new Error('Chosen column must be of Macromolecule semantic type');

        if (state.mode === 'kalign')
          return doKalign(col, table);

        return doEngineMsa(col, table);
      }

      async function doKalign(col: DG.Column<string>, df: DG.DataFrame): Promise<DG.Column<string>> {
        const unusedName = df.columns.getUnusedName(`msa(${col.name})`);
        const sh = seqHelper.getSeqHandler(col);
        const fastaCol = sh.isFasta() ? col : sh.convert(NOTATION.FASTA);
        return runKalign(
          df, fastaCol, false, unusedName, clustersColInput.value,
          kalignGapOpen.value ?? undefined, kalignGapExtend.value ?? undefined,
          kalignTerminalGap.value ?? undefined, onlySelectedInput.value,
        );
      }

      async function doEngineMsa(col: DG.Column<string>, df: DG.DataFrame): Promise<DG.Column<string>> {
        if (!state.currentFunc || !state.currentFuncCall)
          throw new Error('No MSA engine selected');

        // Convert to HELM if needed - prefer notation provider's toHelm if available
        const sh = seqHelper.getSeqHandler(col);
        let srcCol: DG.Column<string>;
        if (sh.isHelm())
          srcCol = col;
        else if (sh.isSeparator() && sh.alphabet === ALPHABET.UN)
          srcCol = sh.convert(NOTATION.HELM);
        else
          srcCol = sh.convert(NOTATION.HELM);

        const func = state.currentFunc;
        const firstParamName = func.inputs[0].name;

        // Read config params from the editor (all params except the first column param)
        const configParams: Record<string, any> = {};
        for (let i = 1; i < func.inputs.length; i++) {
          const name = func.inputs[i].name;
          configParams[name] = state.currentFuncCall.inputs[name];
        }

        const helmResultCol = await runEngineWithClustering(
          func, firstParamName, configParams, srcCol,
          clustersColInput.value, onlySelectedInput.value, df,
        );

        // If column has a notation provider with fromHelm, convert result back to original notation
        const np = col.temp?.[SeqTemps.notationProvider];
        const npCons = np ? np.constructor as typeof NotationProviderBase : null;
        if (npCons?.implementsFromHelm) {
          const convertedName = df.columns.getUnusedName(`msa(${col.name})`);
          const convertedCol = DG.Column.string(convertedName, helmResultCol.length);
          convertedCol.init((i) => {
            const helm = helmResultCol.get(i);
            if (!helm) return '';
            try {
              return npCons.convertFromHelm(helm, {});
            } catch {
              return '';
            }
          });
          convertedCol.semType = DG.SEMTYPE.MACROMOLECULE;
          convertedCol.meta.units = NOTATION.CUSTOM;
          convertedCol.setTag(bioTAGS.aligned, 'SEQ.MSA');
          convertedCol.setTag(bioTAGS.alphabet, ALPHABET.UN);

          // Add HELM column too if requested
          if (includeHelmInput.value)
            df.columns.add(helmResultCol);

          return convertedCol;
        }

        return helmResultCol;
      }

      /** Apply engine and params from options (for programmatic/test use). */
      async function applyEngineOptions(): Promise<void> {
        if (!options.engine || state.mode !== 'engine') return;

        const engine = msaEngines.find(
          (f) => f.name === options.engine || f.friendlyName === options.engine,
        );
        if (!engine) return;

        engineInput.value = engine.friendlyName;
        await updateEngineEditor();
        if (options.engineParams && state.currentFuncCall) {
          for (const [key, value] of Object.entries(options.engineParams))
            state.currentFuncCall.inputs[key] = value;
        }
      }

      // --- Dialog ---

      const dlg = ui.dialog('MSA')
        .add(colInput)
        .add(clustersColInput)
        .add(engineInput)
        .add(engineDescDiv)
        .add(engineParamsButton)
        .add(engineParamsDiv)
        .add(includeHelmInput)
        .add(kalignParamsDiv)
        .add(kalignParamsButton)
        .add(kalignVersionDiv)
        .add(onlySelectedInput)
        .onOK(async () => {
          const pi = DG.TaskBarProgressIndicator.create('Performing MSA...');
          try {
            const resultCol = await doAlignment();
            table.columns.add(resultCol);
            await grok.data.detectSemanticTypes(table);
            if (resultCol.meta.units !== NOTATION.HELM)
              resultCol.setTag(bioTAGS.aligned, 'SEQ.MSA');
            resolve(resultCol);
          } catch (err: any) {
            reject(err);
          } finally {
            pi.close();
          }
        });
      const okBtn = dlg.getButton('OK');

      // Initialize: detect mode from initial column
      switchMode('kalign');
      colInput.fireChanged();

      // If column is pre-specified (tests/programmatic), run immediately without dialog
      if (options.col) {
        await onColumnChanged(options.col);
        await applyEngineOptions();

        const pi = DG.TaskBarProgressIndicator.create('Performing MSA...');
        try {
          const resultCol = await doAlignment();
          table.columns.add(resultCol);
          await grok.data.detectSemanticTypes(table);
          if (resultCol.meta.units !== NOTATION.HELM)
            resultCol.setTag(bioTAGS.aligned, 'SEQ.MSA');
          resolve(resultCol);
        } catch (err: any) {
          reject(err);
        } finally {
          pi.close();
        }
        return;
      }

      dlg.show();
    } catch (err: any) {
      reject(err);
    }
  });
}


/** Runs a discovered MSA engine function with per-cluster alignment support.
 * Groups rows by cluster, creates subset columns, calls the engine per cluster,
 * and merges results into a single output column. */
async function runEngineWithClustering(
  func: DG.Func, colParamName: string, configParams: Record<string, any>,
  srcCol: DG.Column<string>, clustersCol: DG.Column | null,
  onlySelected: boolean, table: DG.DataFrame,
): Promise<DG.Column<string>> {
  const rowCount = srcCol.length;

  // Group rows by cluster
  clustersCol ??= DG.Column.string('Clusters', rowCount).init('0');
  if (clustersCol.type !== DG.COLUMN_TYPE.STRING)
    clustersCol = clustersCol.convertTo(DG.TYPE.STRING);

  const categories = clustersCol.categories;
  const data = clustersCol.getRawData();
  const clusterIndexes: number[][] = new Array(categories.length);

  if (onlySelected) {
    const sel = table.selection;
    if (sel.trueCount === 0)
      throw new Error('No selected rows in the table.');
    for (let i = -1; (i = sel.findNext(i, true)) !== -1;)
      (clusterIndexes[data[i]] ??= []).push(i);
  } else {
    for (let i = 0; i < rowCount; i++)
      (clusterIndexes[data[i]] ??= []).push(i);
  }
  checkForSingleSeqClusters(clusterIndexes, categories);

  const unusedName = table.columns.getUnusedName(`msa(${srcCol.name})`);
  const resultValues: string[] = new Array(rowCount).fill('');
  let lastResultCol: DG.Column<string> | null = null;

  for (const rowIds of clusterIndexes) {
    if (!rowIds || rowIds.length === 0) continue;

    // Create a subset column with just this cluster's sequences
    const subsetSeqs = rowIds.map((i) => srcCol.get(i)!);
    const subsetCol = DG.Column.fromStrings('seq', subsetSeqs);
    copyColumnMetadata(srcCol, subsetCol);
    DG.DataFrame.fromColumns([subsetCol]); // attach to a DataFrame for column operations

    // Call the engine function with the subset
    const call = func.prepare({[colParamName]: subsetCol, ...configParams});
    await call.call();
    const clusterResult = call.getOutputParamValue() as DG.Column<string>;
    lastResultCol = clusterResult;

    // Map cluster results back to original row positions
    for (let i = 0; i < rowIds.length; i++)
      resultValues[rowIds[i]] = clusterResult.get(i) ?? '';
  }

  // Build final column with metadata from the engine's output
  const finalCol = DG.Column.fromStrings(unusedName, resultValues);
  if (lastResultCol) {
    finalCol.meta.units = lastResultCol.meta.units;
    finalCol.semType = lastResultCol.semType;
    for (const tag of [bioTAGS.alphabet, bioTAGS.separator, bioTAGS.alphabetIsMultichar]) {
      const val = lastResultCol.getTag(tag);
      if (val) finalCol.setTag(tag, val);
    }
  }

  return finalCol;
}


function copyColumnMetadata(src: DG.Column, dst: DG.Column): void {
  dst.semType = src.semType;
  dst.meta.units = src.meta.units;
  for (const tag of [bioTAGS.alphabet, bioTAGS.separator, bioTAGS.alphabetIsMultichar]) {
    const val = src.getTag(tag);
    if (val) dst.setTag(tag, val);
  }
}
