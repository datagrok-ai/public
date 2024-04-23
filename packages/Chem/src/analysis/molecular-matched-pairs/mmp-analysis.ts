import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {getRdKitModule}
  from '../../utils/chem-common-rdkit';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ILineSeries, MouseOverLineEvent, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {debounceTime} from 'rxjs/operators';
import {getSigFigs} from '../../utils/chem-common';
import {getMmpFrags, getMmpRules, MmpRules} from './mmp-fragments';
import {getMmpActivityPairsAndTransforms} from './mmp-pairs-transforms';
import {getMmpTrellisPlot} from './mmp-frag-vs-frag';
import {fillPairInfo, getMmpScatterPlot, runMmpChemSpace} from './mmp-cliffs';
import {getInverseSubstructuresAndAlign, PaletteCodes, getPalette} from './mmp-mol-rendering';
import {MMP_COLNAME_FROM, MMP_COLNAME_TO, MMP_COL_PAIRNUM,
  MMP_VIEW_NAME, MMP_TAB_TRANSFORMATIONS, MMP_TAB_FRAGMENTS,
  MMP_TAB_CLIFFS, MMP_TAB_GENERATION, MMP_STRUCT_DIFF_FROM_NAME, MMP_STRUCT_DIFF_TO_NAME} from './mmp-constants';
import {drawMoleculeLabels} from '../../rendering/molecule-label';
import {getGenerations} from './mmp-generations';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import { getEmbeddingColsNames } from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';

export class MmpAnalysis {
  parentTable: DG.DataFrame;
  parentCol: DG.Column;
  mmpRules: MmpRules;
  mmpView: DG.View;
  enableFilters: boolean = true;
  colorPalette: PaletteCodes;
  //transformations tab objects
  allPairsGrid: DG.Grid;
  transformationsMask: DG.BitSet;
  casesGrid: DG.Grid;
  pairsMask: DG.BitSet;
  //cliffs tab objects
  diffs: Array<Float32Array>;
  linesIdxs: Uint32Array;
  cutoffMasks: Array<DG.BitSet>;
  totalCutoffMask: DG.BitSet;
  linesMask: BitArray;
  linesRenderer: ScatterPlotLinesRenderer | null = null;
  lines: ILineSeries;
  linesActivityCorrespondance: Uint32Array;
  //generations tab objects
  generationsGrid: DG.Grid;
  //rdkit
  rdkitModule: RDModule;
  currentTab = '';
  lastSelectedPair: number | null = null;
  propPanelViewer: FormsViewer;

  constructor(table: DG.DataFrame, molecules: DG.Column, palette: PaletteCodes,
    rules: MmpRules, diffs: Array<Float32Array>,
    linesIdxs: Uint32Array, allPairsGrid: DG.Grid, casesGrid: DG.Grid, generationsGrid: DG.Grid,
    tp: DG.Viewer, sp: DG.Viewer,
    sliderInputs: DG.InputBase[], sliderInputValueDivs: HTMLDivElement[], colorInputs: DG.InputBase[],
    linesEditor: ScatterPlotLinesRenderer, lines: ILineSeries, linesActivityCorrespondance: Uint32Array,
    rdkitModule: RDModule) {
    this.rdkitModule = rdkitModule;

    this.parentTable = table;
    this.parentCol = molecules;
    this.colorPalette = palette;
    this.mmpRules = rules;

    this.diffs = diffs;
    this.allPairsGrid = allPairsGrid;
    this.casesGrid = casesGrid;
    this.generationsGrid = generationsGrid;
    this.lines = lines;
    this.linesActivityCorrespondance = linesActivityCorrespondance;
    this.linesIdxs = linesIdxs;

    //transformations tab
    this.transformationsMask = DG.BitSet.create(this.allPairsGrid.dataFrame.rowCount);
    this.transformationsMask.setAll(true);
    this.parentTable.onCurrentRowChanged.pipe(debounceTime(1000)).subscribe(() => {
      if (this.parentTable.currentRowIdx !== -1) {
        this.refilterAllPairs(true);
        this.refreshPair(this.rdkitModule);
      }
    });

    this.refilterAllPairs(true);

    this.allPairsGrid.table.onCurrentRowChanged.subscribe(() => {
      if (this.allPairsGrid.table.currentRowIdx !== -1)
        this.refreshPair(this.rdkitModule);
    });

    this.mmpView = DG.View.create();
    this.mmpView.name = MMP_VIEW_NAME;
    this.mmpView.box = true;

    this.pairsMask = DG.BitSet.create(this.casesGrid.dataFrame.rowCount);
    this.pairsMask.setAll(false);
    this.refreshPair(this.rdkitModule);

    //Cliffs tab
    const roots: any[] = new Array<any>(sliderInputs.length);
    for (let i = 0; i < sliderInputs.length; i ++) {
      sliderInputs[i].onChanged(() => {
        sliderInputValueDivs[i].innerText = sliderInputs[i].value === 0 ? '0' :
          getSigFigs(sliderInputs[i].value, 4).toString();
        this.refilterCliffs(sliderInputs.map((si) => si.value), true);
      });

      ui.tooltip.bind(sliderInputs[i].captionLabel, 'Select the cutoff by activity difference');
      ui.tooltip.bind(sliderInputs[i].input, 'Activity value cutoff');

      colorInputs[i].value = this.colorPalette.hex[i];
      colorInputs[i].onChanged(() => {
        const progressRendering = DG.TaskBarProgressIndicator.create(`Changing colors...`);

        //refresh lines
        for (let i =0; i < colorInputs.length; i++) {
          this.colorPalette.hex[i] = colorInputs[i].value;
          this.colorPalette.numerical[i] = DG.Color.fromHtml(colorInputs[i].value);
          this.colorPalette.rgb[i] = DG.Color.toRgb(this.colorPalette.numerical[i]);
          this.colorPalette.rgbCut[i] = this.colorPalette.rgb[i].replace('rgb(', '').replace(')', '');
        }

        const colors = this.lines.colors;
        for (let i = 0; i < colors!.length; i++)
          colors![i] = this.colorPalette.rgbCut[linesActivityCorrespondance[i]];

        const lines: ILineSeries = {
          from: this.lines.from,
          to: this.lines.to,
          drawArrows: true,
          opacity: 0.5,
          colors: colors,
          arrowSize: 10,
        };

        this.lines = lines;
        this.linesRenderer!.updateLines(lines);

        //refresh trellis plot
        const schemes = new Array<any>(colorInputs.length);
        for (let i = 0; i < colorInputs.length; i++)
          schemes[i] = [this.colorPalette.numerical[i]];

        tp.setOptions({'colorShemes': schemes});
        // tp. ;
        progressRendering.close();
      });

      roots[i] = ui.divV([sliderInputs[i].root, colorInputs[i]]);
    }

    const sliders = ui.divH(roots, 'css-flex-wrap');
    sp.root.style.width = '100%';

    const cliffs1 = ui.divV([sliders, ui.box(sp.root, {style: {maxHeight: '100px', paddingRight: '6px'}})],
      'css-flex-wrap');

    const cliffs = ui.box(cliffs1);

    this.cutoffMasks = new Array<DG.BitSet>(sliderInputs.length);
    this.totalCutoffMask = DG.BitSet.create(this.parentTable.rowCount);
    this.totalCutoffMask.setAll(true);
    this.linesMask = new BitArray(linesIdxs.length);
    this.linesMask.setAll(true);

    for (let i = 0; i < sliderInputs.length; i++) {
      this.cutoffMasks[i] = DG.BitSet.create(this.parentTable.rowCount);
      this.cutoffMasks[i].setAll(true);
    }

    this.linesRenderer = linesEditor;

    this.refilterCliffs(sliderInputs.map((si) => si.value), false);

    //generated tab


    //tabs
    const tabs = ui.tabControl(null, false);

    tabs.addPane(MMP_TAB_TRANSFORMATIONS, () => {

      const addToWorkspaceButton = (table: DG.DataFrame, name: string) => {
        const button = ui.iconFA('arrow-square-down', () => {
          const clonedTable = table.clone();
          clonedTable.name = name;
          grok.shell.addTableView(clonedTable);
        }, 'Add table to workspace');
        button.classList.add('chem-mmpa-add-to-workspace-button');
        return button;
      }

      const createGridDiv = (name: string, grid: DG.Grid) => {
        const header = ui.h1(name, 'chem-mmpa-transformation-tab-header');
        grid.root.prepend(header);
        return ui.splitV([
          ui.box(ui.divH([header, addToWorkspaceButton(grid.dataFrame, name)]), { style: { maxHeight: '30px' } }),
          grid.root
        ]);
      }

      return ui.splitV([
        createGridDiv('Fragments', this.allPairsGrid),
        createGridDiv('Pairs', this.casesGrid),
      ], {}, true);
    });
    tabs.addPane(MMP_TAB_FRAGMENTS, () => {
      return tp.root;
    });
    tabs.addPane(MMP_TAB_CLIFFS, () => {
      return cliffs;
    });
    tabs.addPane(MMP_TAB_GENERATION, () => {
      return this.generationsGrid.root;
    });

    tabs.onTabChanged.subscribe(() => {
      this.currentTab = tabs.currentPane.name;
      if (tabs.currentPane.name == MMP_TAB_TRANSFORMATIONS) {
        this.enableFilters = true;
        this.refilterAllPairs(false);
      } else if (tabs.currentPane.name == MMP_TAB_FRAGMENTS) {
        this.refreshFilterAllPairs();
        this.enableFilters = false;
      } else if (tabs.currentPane.name == MMP_TAB_CLIFFS) {
        this.refilterCliffs(sliderInputs.map((si) => si.value), false);
        if (this.lastSelectedPair) {
          grok.shell.o = fillPairInfo(this.lastSelectedPair, this.linesIdxs,
            this.linesActivityCorrespondance[this.lastSelectedPair],
            this.casesGrid.dataFrame, this.diffs, this.parentTable, this.rdkitModule);
        }
      }
    });

    const decript1 = 'Shows all fragmental substitutions for a given molecule';
    const decript2 = 'Analysis of fragments versus explored value';
    const decript3 = 'Cliffs analysis';

    tabs.getPane(MMP_TAB_TRANSFORMATIONS).header.onmouseover =
      (ev): void => ui.tooltip.show(decript1, ev.clientX, ev.clientY + 5);
    tabs.getPane(MMP_TAB_FRAGMENTS).header.onmouseover =
      (ev): void => ui.tooltip.show(decript2, ev.clientX, ev.clientY + 5);
    tabs.getPane(MMP_TAB_CLIFFS).header.onmouseover =
      (ev): void => ui.tooltip.show(decript3, ev.clientX, ev.clientY + 5);

    this.linesRenderer.lineClicked.subscribe((event: MouseOverLineEvent) => {
      this.linesRenderer!.currentLineId = event.id;
      if (event.id !== -1) {
        grok.shell.o = fillPairInfo(event.id, linesIdxs, linesActivityCorrespondance[event.id],
          casesGrid.dataFrame, diffs, table, rdkitModule, this.propPanelViewer);
        this.lastSelectedPair = event.id;
      }
    });

    this.mmpView.append(tabs);

    const propertiesColumnsNames = this.parentTable.columns.names()
      .filter((name) => name !== this.parentCol.name && !name.startsWith('~'));
    this.propPanelViewer = new FormsViewer();
    this.propPanelViewer.dataframe = this.parentTable;
    this.propPanelViewer.columns = propertiesColumnsNames;
  }

  static async init(table: DG.DataFrame, molecules: DG.Column, activities: DG.ColumnList) {
    //rdkit module
    const module = getRdKitModule();

    //initial calculations
    const frags = await getMmpFrags(molecules);
    const [mmpRules, allCasesNumber] = getMmpRules(frags);
    const palette = getPalette(activities.length);

    //Transformations tab
    const {maxActs, diffs, activityMeanNames, linesIdxs, allPairsGrid, casesGrid, lines, linesActivityCorrespondance} =
      await getMmpActivityPairsAndTransforms(molecules, activities, mmpRules, allCasesNumber, palette);

    //Fragments tab
    const tp = getMmpTrellisPlot(allPairsGrid, activityMeanNames, palette);

    const embedColsNames = getEmbeddingColsNames(table).map((it) => `~${it}`);
    //Cliffs tab
    const [sp, sliderInputs, sliderInputValueDivs, colorInputs] = getMmpScatterPlot(table, activities, maxActs, embedColsNames);
    drawMoleculeLabels(table, molecules, sp as DG.ScatterPlotViewer, 20, 7, 100, 110);

    //running internal chemspace
    const linesEditor = runMmpChemSpace(table, molecules, sp, lines, linesIdxs, linesActivityCorrespondance,
      casesGrid.dataFrame, diffs, module, embedColsNames);

    const generationsGrid: DG.Grid =
      getGenerations(molecules, frags, allPairsGrid, activityMeanNames, activities, module);

    return new MmpAnalysis(table, molecules, palette, mmpRules, diffs, linesIdxs,
      allPairsGrid, casesGrid, generationsGrid,
      tp, sp, sliderInputs, sliderInputValueDivs, colorInputs, linesEditor, lines, linesActivityCorrespondance, module);
  }

  findSpecificRule(diffFromSubstrCol: DG.Column): [idxPairs: number, cases: number[]] {
    const idxParent = this.parentTable.currentRowIdx;
    let idxPairs = -1;
    const cases: number[] = [];
    const idx = this.allPairsGrid.table.currentRowIdx;
    if (idx !== -1) {
      const ruleSmi1 = this.allPairsGrid.table.getCol('From').get(idx);
      const ruleSmi2 = this.allPairsGrid.table.getCol('To').get(idx);
      const ruleSmiNum1 = this.mmpRules.smilesFrags.indexOf(ruleSmi1);
      const ruleSmiNum2 = this.mmpRules.smilesFrags.indexOf(ruleSmi2);

      let counter = 0;

      for (let i = 0; i < this.mmpRules.rules.length; i++) {
        const first = this.mmpRules.rules[i].smilesRule1;
        const second = this.mmpRules.rules[i].smilesRule2;
        for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
          if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
            this.pairsMask.set(counter, true, false);
            if (diffFromSubstrCol.get(counter) === null)
              cases.push(counter);
            if (this.mmpRules.rules[i].pairs[j].firstStructure == idxParent)
              idxPairs = counter;
          }
          counter++;
        }
      }
    }
    return [idxPairs, cases];
  }

  async recoverHighlights(
    cases: number[], diffFrom : DG.Column, diffTo: DG.Column,
    diffFromSubstrCol: DG.Column, diffToSubstrCol: DG.Column, rdkitModule: RDModule) : Promise<void> {
    const pairsFrom = Array<string>(cases.length);
    const pairsTo = Array<string>(cases.length);
    for (let i = 0; i < cases.length; i++) {
      pairsFrom[i] = diffFrom.get(cases[i]);
      pairsTo[i] = diffTo.get(cases[i]);
    }
    const {inverse1, inverse2, fromAligned, toAligned} =
      await getInverseSubstructuresAndAlign(pairsFrom, pairsTo, rdkitModule);
    for (let i = 0; i < cases.length; i++) {
      diffFrom.set(cases[i], fromAligned[i]);
      diffTo.set(cases[i], toAligned[i]);
      diffFromSubstrCol.set(cases[i], inverse1[i]);
      diffToSubstrCol.set(cases[i], inverse2[i]);
    }
  }

  /**
  * Prepares all the entities to show for selected pair.
  * @param {RDModule} rdkitModule RDkit module instance
  */
  async refreshPair(rdkitModule: RDModule) : Promise<void> {
    const progressBarPairs = DG.TaskBarProgressIndicator.create(`Refreshing pairs...`);

    this.pairsMask.setAll(false);
    const diffFromSubstrCol = this.casesGrid.dataFrame.getCol(MMP_STRUCT_DIFF_FROM_NAME);
    const diffToSubstrCol = this.casesGrid.dataFrame.getCol(MMP_STRUCT_DIFF_TO_NAME);
    const diffFrom = this.casesGrid.dataFrame.getCol(MMP_COLNAME_FROM);
    const diffTo = this.casesGrid.dataFrame.getCol(MMP_COLNAME_TO);

    const [idxPairs, cases] = this.findSpecificRule(diffFromSubstrCol);
    await this.recoverHighlights(cases, diffFrom, diffTo, diffFromSubstrCol, diffToSubstrCol, rdkitModule);

    this.casesGrid.dataFrame.filter.copyFrom(this.pairsMask);

    this.casesGrid.setOptions({
      pinnedRowColumnNames: [],
      pinnedRowValues: [],
    });
    if (idxPairs >= 0) {
      this.casesGrid.setOptions({
        pinnedRowValues: [idxPairs.toString()],
        pinnedRowColumnNames: [MMP_COL_PAIRNUM],
      });
    }

    progressBarPairs.close();
  }

  refreshFilterAllPairs() {
    const consistsBitSet: DG.BitSet = DG.BitSet.create(this.allPairsGrid.dataFrame.rowCount);
    consistsBitSet.setAll(true);
    this.allPairsGrid.dataFrame.filter.copyFrom(consistsBitSet);
  }

  /**
  * Gets visible all pairs according to selected molecule in the table.
  * @param {boolean} rowChanged flag if row was changed
  */
  refilterAllPairs(rowChanged: boolean) {
    let idxTrue = -1;
    if (rowChanged) {
      const idx = this.parentTable.currentRowIdx;
      this.transformationsMask.setAll(false);

      for (let i = 0; i < this.mmpRules.rules.length; i++) {
        for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
          const fs = this.mmpRules.rules[i].pairs[j].firstStructure;
          if (idx == fs) {
            if (idxTrue == -1)
              idxTrue = i;
            this.transformationsMask.set(i, true, false);
            break;
          }
        }
      }
    }

    if (this.enableFilters) {
      this.allPairsGrid.dataFrame.filter.copyFrom(this.transformationsMask);
      this.allPairsGrid.invalidate();
      if (rowChanged)
        this.allPairsGrid.table.currentRowIdx = idxTrue;
    }
  }

  /**
  * Makes the dataset filtering according to user specified cutoofs for each activity.
  * @param {number[]} cutoffs for each activity the corresponding cuttof is specified by user
  * @param {boolean} refilter flag if calculation is required before filtering
  */
  refilterCliffs(cutoffs: number[], refilter: boolean) {
    if (refilter) {
      for (let i = 0; i < this.cutoffMasks.length; i ++)
        this.cutoffMasks[i].setAll(false);

      this.totalCutoffMask.setAll(false);
      this.linesMask.setAll(false);

      //setting activity associated masks
      for (let i = 0; i < this.lines.from.length; i++) {
        const activityNumber = this.linesActivityCorrespondance[i];
        const line = this.linesIdxs[i];
        if (this.diffs[activityNumber][line] >= cutoffs[activityNumber]) {
          this.cutoffMasks[activityNumber].set(this.lines.from[i], true, false);
          this.cutoffMasks[activityNumber].set(this.lines.to[i], true, false);
          this.linesMask.setBit(i, true, false);
        }
      }

      //setting resulting masks
      //TODO: refine
      for (let j = 0; j < this.totalCutoffMask.length; j++) {
        let res = false;
        for (let i = 0; i < this.cutoffMasks.length; i++)
          res = res || this.cutoffMasks[i].get(j);
        this.totalCutoffMask.set(j, res, false);
      }
    }

    this.parentTable.filter.copyFrom(this.totalCutoffMask);
    this.linesRenderer!.linesVisibility = this.linesMask;
  }
}
