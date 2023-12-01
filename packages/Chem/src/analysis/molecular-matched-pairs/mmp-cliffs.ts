import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ISubstruct} from '../../rendering/rdkit-cell-renderer';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {convertMolNotation} from '../../package';
import {drawMoleculeToCanvas} from '../../utils/chem-common-rdkit';
import {getSigFigs} from '../../utils/chem-common';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {ILineSeries, MouseOverLineEvent, ScatterPlotCurrentLineStyle, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {chemSpace} from '../chem-space';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {debounceTime} from 'rxjs/operators';
import {getInverseSubstructuresAndAlign} from './mmp-mol-rendering';
import {MMP_COLNAME_FROM, MMP_COLNAME_TO, MMP_COL_PAIRNUM_FROM, MMP_COL_PAIRNUM_TO, MMP_COLNAME_CHEMSPACE_X,
  MMP_COLNAME_CHEMSPACE_Y, MMP_STRUCT_DIFF_FROM_NAME, MMP_STRUCT_DIFF_TO_NAME} from './mmp-constants';

export let lastSelectedPair: number | null = null;

export function getMmpScatterPlot(table: DG.DataFrame, activities: DG.ColumnList, maxActs: number[]) :
[sp: DG.Viewer, sliderInputs: DG.InputBase[], sliderInputValueDivs: HTMLDivElement[], colorInputs: DG.InputBase[]] {
  const colX = DG.Column.float('~X', table.rowCount);
  const colY = DG.Column.float('~Y', table.rowCount);
  table.columns.add(colX);
  table.columns.add(colY);
  const sp = DG.Viewer.scatterPlot(table, {
    x: '~X',
    y: '~Y',
    zoomAndFilter: 'no action',
    //color: activities.name,
    showXSelector: false,
    showYSelector: false,
    markerDefaultSize: 7,
    markerType: 'circle border',
  });

  const sliderInputs = new Array<DG.InputBase>(maxActs.length);
  const sliderInputValueDivs = new Array<HTMLDivElement>(maxActs.length);
  const colorInputs = new Array<DG.InputBase>(maxActs.length);

  for (let i = 0; i < maxActs.length; i ++) {
    const sliderInput = ui.sliderInput(activities.byIndex(i).name, 0, 0, maxActs[i]);
    sliderInput.root.style.marginLeft = '6px';
    const sliderInputValueDiv = ui.divText(sliderInput.stringValue, 'ui-input-description');
    sliderInput.addOptions(sliderInputValueDiv);
    sliderInputs[i] = sliderInput;
    sliderInputValueDivs[i] = sliderInputValueDiv;
    const colorInput = ui.colorInput('', '#FF0000');
    colorInput.root.style.marginLeft = '6px';
    colorInputs[i] = colorInput;
  }

  return [sp, sliderInputs, sliderInputValueDivs, colorInputs];
}

function drawMolPair(molecules: string[], substr: (ISubstruct | null)[], div: HTMLDivElement, tooltip?: boolean) {
  ui.empty(div);
  const hosts = ui.divH([]);
  if (!tooltip)
    hosts.classList.add('chem-mmpa-context-pane-mol-div');
  const canvasWidth = 170;
  const canvasHeight = 75;
  for (let i = 0; i < 2; i++) {
    const imageHost = ui.canvas(canvasWidth, canvasHeight);
    drawMoleculeToCanvas(0, 0, canvasWidth, canvasHeight, imageHost, molecules[i], '', undefined, substr[i]);
    hosts.append(imageHost);
  }
  div.append(hosts);
};

export function moleculesPairInfo(line: number, linesIdxs: Uint32Array, activityNum: number, pairsDf: DG.DataFrame,
  diffs: Array<Float32Array>, parentTable: DG.DataFrame, molColName: string, rdkitModule: RDModule, tooltip?: boolean):
  HTMLDivElement {
  const div = ui.divV([], {style: {width: '100%'}});
  const moleculesDiv = ui.divH([]);
  div.append(moleculesDiv);
  const pairIdx = linesIdxs[line];
  const subsrtFrom = pairsDf.get(MMP_STRUCT_DIFF_FROM_NAME, pairIdx);
  const subsrtTo = pairsDf.get(MMP_STRUCT_DIFF_TO_NAME, pairIdx);
  const moleculeFrom = pairsDf.get(MMP_COLNAME_FROM, pairIdx);
  const moleculeTo = pairsDf.get(MMP_COLNAME_TO, pairIdx);
  if (!tooltip) {
    const fromIdx = pairsDf.get(MMP_COL_PAIRNUM_FROM, pairIdx);
    const toIdx = pairsDf.get(MMP_COL_PAIRNUM_TO, pairIdx);
    const props = getMoleculesPropertiesDiv([fromIdx, toIdx], parentTable, molColName);
    div.append(props);
  } else {
    const diff = ui.tableFromMap({'Diff': getSigFigs(diffs[activityNum][pairIdx], 4)});
    diff.style.maxWidth = '150px';
    div.append(diff);
  }
  if (subsrtFrom || subsrtTo)
    drawMolPair([moleculeFrom, moleculeTo], [subsrtFrom, subsrtTo], moleculesDiv, tooltip);
  else {
    moleculesDiv.append(ui.divText(`Loading...`));
    getInverseSubstructuresAndAlign([moleculeFrom], [moleculeTo], rdkitModule).then((res) => {
      const {inverse1, inverse2, fromAligned, toAligned} = res;
      pairsDf.set(MMP_STRUCT_DIFF_FROM_NAME, pairIdx, inverse1[0]);
      pairsDf.set(MMP_STRUCT_DIFF_TO_NAME, pairIdx, inverse2[0]);
      pairsDf.set(MMP_COLNAME_FROM, pairIdx, fromAligned[0]);
      pairsDf.set(MMP_COLNAME_TO, pairIdx, toAligned[0]);
      drawMolPair([fromAligned[0], toAligned[0]], [inverse1[0], inverse2[0]], moleculesDiv, tooltip);
    });
  }
  return div;
};

function getMoleculesPropertiesDiv(idxs: number[], parentTable: DG.DataFrame, molColName: string): HTMLElement {
  const propertiesColumnsNames = parentTable.columns.names()
    .filter((name) => name !== molColName && name !== MMP_COLNAME_CHEMSPACE_X &&
      name !== MMP_COLNAME_CHEMSPACE_Y && !name.startsWith('~'));
  const formsViewer = new FormsViewer();
  //@ts-ignore
  formsViewer.dataframe = parentTable;
  formsViewer.columns = propertiesColumnsNames;
  formsViewer.fixedRowNumbers = idxs;
  formsViewer.root.classList.add('chem-mmpa-forms-viewer');

  return ui.div(formsViewer.root, {style: {height: '100%'}});
}

export function runMmpChemSpace(table: DG.DataFrame, molecules: DG.Column, sp: DG.Viewer, lines: ILineSeries,
  linesIdxs: Uint32Array, linesActivityCorrespondance: Uint32Array, pairsDf: DG.DataFrame, diffs: Array<Float32Array>,
  rdkitModule: RDModule): ScatterPlotLinesRenderer {
  const chemSpaceParams = {
    seqCol: molecules,
    methodName: DimReductionMethods.UMAP,
    similarityMetric: BitArrayMetricsNames.Tanimoto as BitArrayMetrics,
    embedAxesNames: [MMP_COLNAME_CHEMSPACE_X, MMP_COLNAME_CHEMSPACE_Y],
    options: {},
  };

  //@ts-ignore
  const spEditor = new ScatterPlotLinesRenderer(sp as DG.ScatterPlotViewer,
    MMP_COLNAME_CHEMSPACE_X, MMP_COLNAME_CHEMSPACE_Y,
    lines, ScatterPlotCurrentLineStyle.bold);


  spEditor.lineClicked.subscribe((event: MouseOverLineEvent) => {
    spEditor.currentLineId = event.id;
    if (event.id !== -1) {
      grok.shell.o = moleculesPairInfo(event.id, linesIdxs, linesActivityCorrespondance[event.id],
        pairsDf, diffs, table, molecules.name, rdkitModule);
      lastSelectedPair = event.id;
    }
  });

  spEditor.lineHover.pipe(debounceTime(500)).subscribe((event: MouseOverLineEvent) => {
    ui.tooltip.show(
      moleculesPairInfo(event.id, linesIdxs, linesActivityCorrespondance[event.id],
        pairsDf, diffs, table, molecules.name, rdkitModule, true),
      event.x, event.y);
  });

  const progressBarSpace = DG.TaskBarProgressIndicator.create(`Running Chemical space...`);
  chemSpace(chemSpaceParams).then((res) => {
    const embeddings = res.coordinates;
    for (const col of embeddings)
      table.columns.replace(col.name, col);
    progressBarSpace.close();
  });

  return spEditor;
}
