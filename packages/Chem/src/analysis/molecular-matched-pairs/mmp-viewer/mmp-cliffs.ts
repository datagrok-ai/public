import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {drawMoleculeToCanvas} from '../../../utils/chem-common-rdkit';
import {getSigFigs} from '../../../utils/chem-common';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {ILineSeries, MouseOverLineEvent, ScatterPlotCurrentLineStyle, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {debounceTime} from 'rxjs/operators';
import {getInverseSubstructuresAndAlign} from './mmp-mol-rendering';
import {MMP_CONTEXT_PANE_CLASS, MMP_NAMES} from './mmp-constants';
import $ from 'cash-dom';
import {MMPA} from '../mmp-analysis/mmpa';
import {ISequenceSpaceParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

export function getMmpScatterPlot(
  parentTable: DG.DataFrame, axesColsNames: string[], labelsColName: string) : DG.ScatterPlotViewer {
  parentTable.columns.addNewFloat(axesColsNames[0]);
  parentTable.columns.addNewFloat(axesColsNames[1]);
  const sp = DG.Viewer.scatterPlot(parentTable, {
    x: axesColsNames[0],
    y: axesColsNames[1],
    zoomAndFilter: 'no action',
    //color: activities.name,
    showXSelector: false,
    showXAxis: false,
    showYSelector: false,
    showYAxis: false,
    showColorSelector: false,
    showSizeSelector: false,
    markerDefaultSize: 7,
    showMarkerLabels: 'Auto',
  });
  //temporary fix (to save backward compatibility) since labels
  //option type has been changed from string to array in 1.23 platform version
  const spProps = Object.keys(sp.props);
  if (spProps.includes('labelColumnNames')) { //@ts-ignore
    if (sp.props['labelColumnNames'].constructor.name == 'Array')
      sp.setOptions({labelColumnNames: [labelsColName]});
  }
  if (spProps.includes('useLabelAsMarker')) { //@ts-ignore
    sp.setOptions({useLabelAsMarker: true, labelAsMarkerSize: 50});
  }
  return sp;
}

function drawMolPair(molecules: string[], indexes: number[], substruct: (ISubstruct | null)[], div: HTMLDivElement,
  parentTable: DG.DataFrame, tooltip?: boolean) {
  ui.empty(div);
  const hosts = ui.divH([]);
  if (!tooltip)
    hosts.classList.add('chem-mmpa-context-pane-mol-div');
  const canvasWidth = 170;
  const canvasHeight = 75;
  for (let i = 0; i < 2; i++) {
    const imageHost = ui.canvas(canvasWidth, canvasHeight);
    imageHost.onclick = () => {parentTable.currentRowIdx = indexes[i];};
    drawMoleculeToCanvas(0, 0, canvasWidth, canvasHeight, imageHost, molecules[i], '', undefined, substruct[i]);
    hosts.append(imageHost);
  }
  div.append(hosts);
};

export function fillPairInfo(mmpa: MMPA, line: number, linesIdxs: Uint32Array, activityNum: number,
  pairsDf: DG.DataFrame, diffs: Array<Float32Array>, parentTable: DG.DataFrame,
  rdkitModule: RDModule, propPanelViewer?: FormsViewer):
  HTMLDivElement {
  const div = ui.divV([]);
  $(div).css({'width': '100%', 'height': '100%'});
  const moleculesDiv = ui.divH([]);
  div.append(moleculesDiv);
  const pairIdx = linesIdxs[line];
  const subsrtFrom = pairsDf.get(MMP_NAMES.STRUCT_DIFF_FROM_NAME, pairIdx);
  const subsrtTo = pairsDf.get(MMP_NAMES.STRUCT_DIFF_TO_NAME, pairIdx);
  const moleculeFrom = pairsDf.get(MMP_NAMES.FROM, pairIdx);
  const moleculeTo = pairsDf.get(MMP_NAMES.TO, pairIdx);
  const fromIdx = pairsDf.get(MMP_NAMES.PAIRNUM_FROM, pairIdx);
  const toIdx = pairsDf.get(MMP_NAMES.PAIRNUM_TO, pairIdx);
  const ruleNum = pairsDf.get(MMP_NAMES.RULENUM, pairIdx);
  if (propPanelViewer) {
    const props = getMoleculesPropertiesDiv(propPanelViewer, [fromIdx, toIdx]);
    div.append(props);
  } else {
    const diff = ui.tableFromMap({'Diff': getSigFigs(diffs[activityNum][pairIdx], 4)});
    diff.style.maxWidth = '150px';
    div.append(diff);
    if (subsrtFrom || subsrtTo) {
      drawMolPair([moleculeFrom, moleculeTo], [fromIdx, toIdx],
        [subsrtFrom, subsrtTo], moleculesDiv, parentTable, !propPanelViewer);
    } else {
      moleculesDiv.append(ui.divText(`Loading...`));
      const cores = mmpa.rules.rules[ruleNum].pairs
        .filter((pair) => pair.firstStructure === fromIdx && pair.secondStructure === toIdx);
      if (cores.length) {
        getInverseSubstructuresAndAlign([cores[0].core], [moleculeFrom], [moleculeTo], rdkitModule).then((res) => {
          const {inverse1, inverse2, fromAligned, toAligned} = res;
          pairsDf.set(MMP_NAMES.STRUCT_DIFF_FROM_NAME, pairIdx, inverse1[0]);
          pairsDf.set(MMP_NAMES.STRUCT_DIFF_TO_NAME, pairIdx, inverse2[0]);
          pairsDf.set(MMP_NAMES.FROM, pairIdx, fromAligned[0]);
          pairsDf.set(MMP_NAMES.TO, pairIdx, toAligned[0]);
          drawMolPair([fromAligned[0], toAligned[0]], [fromIdx, toIdx],
            [inverse1[0], inverse2[0]], moleculesDiv, parentTable, !!propPanelViewer);
        });
      }
    }
  }
  div.classList.add(MMP_CONTEXT_PANE_CLASS);
  return div;
};

function getMoleculesPropertiesDiv(propPanelViewer: FormsViewer, idxs: number[]): HTMLElement {
  propPanelViewer.fixedRowNumbers = idxs;
  propPanelViewer.root.classList.add('chem-mmpa-forms-viewer');

  return ui.div(propPanelViewer.root, {style: {height: '100%'}});
}

export function runMmpChemSpace(parentTable: DG.DataFrame, molCol: DG.Column, sp: DG.Viewer, lines: ILineSeries,
  linesIdxs: Uint32Array, linesActivityCorrespondance: Uint32Array, pairsDf: DG.DataFrame, mmpa: MMPA,
  rdkitModule: RDModule, embedColsNames: string[]): [ScatterPlotLinesRenderer, ISequenceSpaceParams] {
  const chemSpaceParams = {
    seqCol: molCol,
    methodName: DimReductionMethods.UMAP,
    similarityMetric: BitArrayMetricsNames.Tanimoto as BitArrayMetrics,
    embedAxesNames: embedColsNames,
    options: {useWebGPU: true},
  };

  const spEditor = new ScatterPlotLinesRenderer(sp as DG.ScatterPlotViewer,
    embedColsNames[0], embedColsNames[1], lines, ScatterPlotCurrentLineStyle.bold);

  spEditor.lineHover.pipe(debounceTime(500)).subscribe((event: MouseOverLineEvent) => {
    ui.tooltip.show(
      fillPairInfo(mmpa, event.id, linesIdxs, linesActivityCorrespondance[event.id],
        pairsDf, mmpa.allCasesBased.diffs, parentTable, rdkitModule),
      event.x, event.y);
  });

  return [spEditor, chemSpaceParams];
}
