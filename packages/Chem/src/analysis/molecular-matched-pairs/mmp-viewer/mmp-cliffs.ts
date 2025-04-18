import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {drawMoleculeToCanvas} from '../../../utils/chem-common-rdkit';
import {getSigFigs} from '../../../utils/chem-common';
import {ILineSeries, MouseOverLineEvent, ScatterPlotCurrentLineStyle, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {debounceTime} from 'rxjs/operators';
import {getInverseSubstructuresAndAlign} from './mmp-mol-rendering';
import {MMP_CONTEXT_PANE_CLASS, MMP_NAMES, MOL_CANVAS_HEIGHT, MOL_CANVAS_WIDTH} from './mmp-constants';
import $ from 'cash-dom';
import {MMPA} from '../mmp-analysis/mmpa';
import {ISequenceSpaceParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

export function getMmpScatterPlot(
  parentTable: DG.DataFrame, axesColsNames: string[], labelsColName: string, colorByActivityName: string) : DG.ScatterPlotViewer {
  parentTable.columns.addNewFloat(axesColsNames[0]);
  parentTable.columns.addNewFloat(axesColsNames[1]);
  const sp = DG.Viewer.scatterPlot(parentTable, {
    x: axesColsNames[0],
    y: axesColsNames[1],
    zoomAndFilter: 'no action',
    color: colorByActivityName,
    showXSelector: false,
    showXAxis: false,
    showYSelector: false,
    showYAxis: false,
    showColorSelector: true,
    showSizeSelector: false,
    markerDefaultSize: 4,
    displayLabels: 'Auto',
    markerType: 'circle',
    jitterSize: 5,
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

function drawMolPairTooltop(molecules: string[], substruct: (ISubstruct | null)[], div: HTMLDivElement) {
  ui.empty(div);
  const hosts = ui.divH([]);
  const molHosts = getMoleculePairImages(molecules, substruct);
  for (const imageHost of molHosts)
    hosts.append(imageHost);
  div.append(hosts);
};

function drawMolPairContextPane(molecules: string[], indexes: number[], substruct: (ISubstruct | null)[],
  divArray: HTMLDivElement[],
  parentTable: DG.DataFrame) {
  divArray.forEach((it) => ui.empty(it));
  const molHosts = getMoleculePairImages(molecules, substruct, parentTable, indexes);
  divArray[0].append(molHosts[0]);
  divArray[1].append(molHosts[1]);
};

export function getMoleculePairImages(molecules: string[], substruct: (ISubstruct | null)[],
  parentTable?: DG.DataFrame, indexes?: number[]): HTMLCanvasElement[] {
  const molHosts = [];
  for (let i = 0; i < 2; i++) {
    const imageHost = ui.canvas(MOL_CANVAS_WIDTH, MOL_CANVAS_HEIGHT);
    if (indexes && parentTable)
      imageHost.onclick = () => {parentTable.currentRowIdx = indexes[i];};
    drawMoleculeToCanvas(0, 0, MOL_CANVAS_WIDTH, MOL_CANVAS_HEIGHT, imageHost,
      molecules[i], '', undefined, substruct[i]);
    molHosts.push(imageHost);
  }
  return molHosts;
}

export function fillPairInfo(mmpa: MMPA, pairIdx: number, pairsDf: DG.DataFrame,
  parentTable: DG.DataFrame, rdkitModule: RDModule, molColName: string, activityNum?: number,
  diffs?: Array<Float32Array>):
  HTMLDivElement {
  const moleculeFrom = pairsDf.get(MMP_NAMES.FROM, pairIdx);
  const moleculeTo = pairsDf.get(MMP_NAMES.TO, pairIdx);
  const fromIdx = pairsDf.get(MMP_NAMES.PAIRNUM_FROM, pairIdx);
  const toIdx = pairsDf.get(MMP_NAMES.PAIRNUM_TO, pairIdx);

  const div = activityNum == undefined ? ui.divH([]) : ui.divV([]);
  let moleculesDiv: HTMLDivElement | null = null;
  let contextPaneMolDivs: HTMLDivElement[] | null = null;
  if (activityNum == undefined) {
    contextPaneMolDivs = [ui.div(ui.divText(`Loading...`)), ui.div(ui.divText(`Loading...`))];
    const headers = ui.divV([], 'chem-mmp-context-pane-headers');
    const mol1Div = ui.divV([]);
    const mol2Div = ui.divV([]);
    const grid = grok.shell.tableView(parentTable.name)?.grid;
    const propertiesColumnsNames = parentTable!.columns.names()
      .filter((name) => !name.startsWith('~'));
    for (const colName of propertiesColumnsNames) {
      if (colName === molColName) {
        const molColHeader = ui.div(ui.divText(colName, 'chem-mmp-context-pane-item'),
          {style: {height: `${MOL_CANVAS_HEIGHT}px`}});
        molColHeader.classList.add('chem-mmp-context-pane-mol-header');
        headers.append(molColHeader);
        mol1Div.append(contextPaneMolDivs[0]);
        mol2Div.append(contextPaneMolDivs[1]);
      } else {
        if (grid) {
          headers.append(ui.divText(colName, 'chem-mmp-context-pane-item'));
          setContextPanelValue(colName, fromIdx, parentTable, grid, mol1Div);
          setContextPanelValue(colName, toIdx, parentTable, grid, mol2Div);
        } else
          grok.shell.error('Parent grid not found');
      }
    }
    div.append(headers, mol1Div, mol2Div);
  } else {
    $(div).css({'width': '100%', 'height': '100%'});
    moleculesDiv = ui.divH([]);
    div.append(moleculesDiv);
    const diff = ui.tableFromMap({'Diff': getSigFigs(diffs![activityNum][pairIdx], 4)});
    diff.style.maxWidth = '150px';
    div.append(diff);
    moleculesDiv.append(ui.divText(`Loading...`));
  }
  const core = mmpa.frags.idToName[pairsDf.get(MMP_NAMES.CORE_NUM, pairIdx)].replace('[*:1]', '[H]');
  if (core) {
    getInverseSubstructuresAndAlign([core],
      [moleculeFrom], [moleculeTo], rdkitModule).then((res) => {
      const {inverse1, inverse2, fromAligned, toAligned} = res;

      if (activityNum != undefined && moleculesDiv) //tooptip
        drawMolPairTooltop([fromAligned[0], toAligned[0]], [inverse1[0], inverse2[0]], moleculesDiv);
      else { //context panel
        if (contextPaneMolDivs) {
          drawMolPairContextPane([fromAligned[0], toAligned[0]], [fromIdx, toIdx],
            [inverse1[0], inverse2[0]], contextPaneMolDivs, parentTable);
        }
      }
    });
  }
  div.classList.add(MMP_CONTEXT_PANE_CLASS);
  return div;
};

function setContextPanelValue(colName: string, idx: number, parentTable: DG.DataFrame,
  grid: DG.Grid, div: HTMLDivElement) {
  let val = parentTable.get(colName, idx);
  if (typeof val === 'number')
    val = DG.format(val, grid.col(colName)?.format);
  div.append(ui.div(ui.divText(val, 'chem-mmp-context-pane-item')));
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
      fillPairInfo(mmpa, linesIdxs[event.id], pairsDf, parentTable, rdkitModule, molCol.name,
        linesActivityCorrespondance[event.id], mmpa.allCasesBased.diffs),
      event.x, event.y);
  });

  return [spEditor, chemSpaceParams];
}
