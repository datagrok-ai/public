import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {MMP_NAMES, MMP_CONSTRICTIONS, MMP_ERRORS} from '../mmp-constants';
import {MmpRules, MmpInput} from '../mmp-constants';
import {getInverseSubstructuresAndAlign, PaletteCodes, getPalette} from '../mmp-mol-rendering';
import {getMmpFrags, getMmpRules} from '../mmp-fragments';
import {getMmpActivityPairsAndTransforms} from '../mmp-pairs-transforms';
import {getMmpTrellisPlot} from '../mmp-frag-vs-frag';
import {getMmpScatterPlot, runMmpChemSpace} from '../mmp-cliffs';
import {getGenerations} from '../mmp-generations';


import {drawMoleculeLabels} from '../../../rendering/molecule-label';
import {ILineSeries, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {getEmbeddingColsNames} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import $ from 'cash-dom';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {getMmpFilters, MmpFilters} from '../mmp-filters';
import {Subject} from 'rxjs/internal/Subject';


export class MatchedMolecularPairsViewer extends DG.JsViewer {
  static TYPE: string = 'MMP';

  //properties
  molecules: string | null = null;
  activities: string[] | null = null;
  fragmentCutoff: number | null;


  //saved
  totalData: string;

  totalDataUpdated: boolean = false;

  moleculesCol: DG.Column | null = null;

  activitiesCols: DG.ColumnList | null = null;

  parentTable: DG.DataFrame | null = null;
  parentCol: DG.Column| null = null;
  mmpRules: MmpRules | null = null;
  mmpView: DG.View | null = null;
  enableFilters: boolean = true;
  colorPalette: PaletteCodes | null = null;
  //transformations tab objects
  transFragmentsGrid: DG.Grid | null = null;
  transFragmentsMask: DG.BitSet | null = null;
  transPairsGrid: DG.Grid | null = null;
  transPairsMask: DG.BitSet | null = null;
  //cliffs tab objects
  diffs: Array<Float32Array> | null = null;
  linesIdxs: Uint32Array | null = null;
  cutoffMasks: Array<DG.BitSet> | null = null;
  totalCutoffMask: DG.BitSet | null = null;
  linesMask: BitArray | null = null;
  linesRenderer: ScatterPlotLinesRenderer | null = null;
  lines: ILineSeries | null = null;
  linesActivityCorrespondance: Uint32Array | null = null;
  //generations tab objects
  generationsGrid: DG.Grid | null = null;
  //rdkit
  rdkitModule: RDModule | null = null;
  currentTab = '';
  lastSelectedPair: number | null = null;
  propPanelViewer: FormsViewer | null = null;
  mmpFilters: MmpFilters | null = null;
  calculatedOnGPU: boolean | null = null;

  constructor() {
    super();
    DG.debounce(this.onPropertyChangedObs, 50).subscribe(this.onPropertyChangedDebounced.bind(this));
    //properties
    this.molecules = this.string('molecules');
    this.activities = this.stringList('activities');
    this.fragmentCutoff = this.float('fragmentCutoff');

    this.totalData = this.string('totalData', 'null', {userEditable: false});
  }

  onTableAttached() {

  }

  onPropertyChangedDebounced() {
    if (this.totalDataUpdated) {
      this.moleculesCol = this.dataFrame.col(this.molecules!);
      this.activitiesCols = DG.DataFrame.fromColumns(this.dataFrame.columns.byNames(this.activities!)).columns;

      //const aa = JSON.parse(this.totalData);

      this.render();
      this.totalDataUpdated = false;
      return;
    }

    this.moleculesCol = this.dataFrame.col(this.molecules!);
    this.activitiesCols = DG.DataFrame.fromColumns(this.dataFrame.columns.byNames(this.activities!)).columns;
    if (this.molecules && this.activities && this.fragmentCutoff) {
      this.render();
      return;
    }
  }

  onPropertyChangedObs : Subject<DG.Property | null> = new Subject<DG.Property | null>();

  onPropertyChanged(property: DG.Property | null): void {
    console.log(property!.name);
    super.onPropertyChanged(property);
    if (property?.name === 'totalData')
      this.totalDataUpdated = true;

    this.onPropertyChangedObs.next(property);
  }

  async render() {
    $(this.root).empty();
    if (this.dataFrame) {
      const loader = ui.div(ui.loader());
      loader.classList.add('mmpa-loader');
      this.root.appendChild(loader);
      const progressMMP = DG.TaskBarProgressIndicator.create(`Running MMP analysis...`);

      try {
        await this.runMMP(
          {table: this.dataFrame,
            molecules: this.moleculesCol!,
            activities: this.activitiesCols!,
            fragmentCutoff: this.fragmentCutoff!,
          });
      } catch (e: any) {

      } finally {
        $(this.root).empty();
        if (this.mmpView)
          this.root.appendChild(this.mmpView!.root);
        else
          this.close();
        progressMMP.close();
      }

      // let tableInfo = this.parentTable!.getTableInfo();
      // await grok.dapi.tables.uploadDataFrame(this.parentTable!);
      // await grok.dapi.tables.save(tableInfo);
      // console.log(`${grok.dapi.root}/entities/${tableInfo.id}/token`);
      // const response = await fetch(`${grok.dapi.root}/entities/${tableInfo.id}/token`, {method: "POST"});
      //or without POST
      // const respText = await response.text();
    }
  }

  //setup after calculation
  setupTransformationTab!: ()=> void;
  setupFilters!: (mmpFilters: MmpFilters, linesActivityCorrespondance: Uint32Array, tp: DG.Viewer) => void;
  setupCliffsTab!: (sp: DG.Viewer, mmpFilters: MmpFilters, linesEditor: ScatterPlotLinesRenderer) => HTMLDivElement;
  getTabs!: (tp: DG.Viewer, mmpFilters: MmpFilters, cliffs: HTMLDivElement) => DG.TabControl;
  fillAll!: (mmpInput: MmpInput, palette: PaletteCodes,
    rules: MmpRules, diffs: Array<Float32Array>,
    linesIdxs: Uint32Array, transFragmentsGrid: DG.Grid, transPairsGrid: DG.Grid, generationsGrid: DG.Grid,
    tp: DG.Viewer, sp: DG.Viewer, mmpFilters: MmpFilters,
    linesEditor: ScatterPlotLinesRenderer, lines: ILineSeries, linesActivityCorrespondance: Uint32Array,
    rdkitModule: RDModule, gpuUsed: boolean) => void;

  async runMMP(mmpInput: MmpInput) {
    //console.profile('MMP');

    const module = getRdKitModule();
    const moleculesArray = mmpInput.molecules.toList();

    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;

    if (!gpu && mmpInput.molecules.length > MMP_CONSTRICTIONS.CPU) {
      grok.shell.error(MMP_ERRORS.FRAGMENTS_CPU);
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    } else if (mmpInput.molecules.length > MMP_CONSTRICTIONS.GPU) {
      grok.shell.error(MMP_ERRORS.FRAGMENTS_GPU);
      throw new Error(MMP_ERRORS.FRAGMENTS_GPU);
    }

    //initial calculations
    const fragsOut = await getMmpFrags(moleculesArray);
    const [mmpRules, allCasesNumber] = await getMmpRules(fragsOut, mmpInput.fragmentCutoff, gpu);

    const palette = getPalette(mmpInput.activities.length);
    //Transformations tab
    const {maxActs, diffs, meanDiffs, activityMeanNames,
      linesIdxs, transFragmentsGrid, transPairsGrid, lines, linesActivityCorrespondance} =
      getMmpActivityPairsAndTransforms(mmpInput, mmpRules, allCasesNumber, palette);

    //Fragments tab
    const tp = getMmpTrellisPlot(transFragmentsGrid, activityMeanNames, palette);

    //Cliffs tab
    const embedColsNames = getEmbeddingColsNames(mmpInput.table).map((it) => `~${it}`);

    const mmpFilters = getMmpFilters(mmpInput, maxActs, transFragmentsGrid.dataFrame.col(MMP_NAMES.PAIRS)!.stats.max);
    console.log(`created mmpa filters`);

    const sp = getMmpScatterPlot(mmpInput, embedColsNames);

    drawMoleculeLabels(mmpInput.table, mmpInput.molecules, sp as DG.ScatterPlotViewer, 20, 7, 100, 110);

    //running internal chemspace
    const linesEditor = runMmpChemSpace(mmpInput, sp, lines, linesIdxs, linesActivityCorrespondance,
      transPairsGrid.dataFrame, diffs, module, embedColsNames);


    const generationsGrid: DG.Grid =
      await getGenerations(mmpInput, moleculesArray, fragsOut, meanDiffs, transFragmentsGrid, activityMeanNames, gpu);


    this.fillAll(mmpInput, palette, mmpRules, diffs, linesIdxs, transFragmentsGrid, transPairsGrid, generationsGrid,
      tp, sp, mmpFilters, linesEditor, lines, linesActivityCorrespondance, module, gpu);

    this.totalData = JSON.stringify({'lucky key': 'luck'});
    //console.profileEnd('MMP');
  }

  findSpecificRule(diffFromSubstrCol: DG.Column): [idxPairs: number, cases: number[]] {
    const idxParent = this.parentTable!.currentRowIdx;
    let idxPairs = -1;
    const cases: number[] = [];
    const idx = this.transFragmentsGrid!.table.currentRowIdx;
    if (idx !== -1) {
      const ruleSmi1 = this.transFragmentsGrid!.table.getCol(MMP_NAMES.FROM).get(idx);
      const ruleSmi2 = this.transFragmentsGrid!.table.getCol(MMP_NAMES.TO).get(idx);
      const ruleSmiNum1 = this.mmpRules!.smilesFrags.indexOf(ruleSmi1);
      const ruleSmiNum2 = this.mmpRules!.smilesFrags.indexOf(ruleSmi2);

      let counter = 0;

      for (let i = 0; i < this.mmpRules!.rules.length; i++) {
        const first = this.mmpRules!.rules[i].smilesRule1;
        const second = this.mmpRules!.rules[i].smilesRule2;
        for (let j = 0; j < this.mmpRules!.rules[i].pairs.length; j++) {
          if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
            this.transPairsMask!.set(counter, true, false);
            if (diffFromSubstrCol.get(counter) === null)
              cases.push(counter);
            if (this.mmpRules!.rules[i].pairs[j].firstStructure == idxParent)
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

  //interaction logics
  refreshPair!: (rdkitModule: RDModule) => Promise<void>;
  refreshFilterAllFragments!: () => void;
  refilterAllFragments!: (rowChanged: boolean) => void;
  refilterCliffs!: (cutoffs: number[], isActiveVar: boolean[], refilter: boolean) => void;
  pinPair!: (idx: number) => void;
}

import(/* webpackMode: "eager" */ './mmp-viewer-setup');
import(/* webpackMode: "eager" */ './mmp-viewer-transformation-related');
import(/* webpackMode: "eager" */ './mmp-viewer-cliffs-related');
