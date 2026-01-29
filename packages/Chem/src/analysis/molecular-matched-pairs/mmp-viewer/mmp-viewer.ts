import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {Subject} from 'rxjs/internal/Subject';
import $ from 'cash-dom';
import {ILineSeries, MouseOverLineEvent, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';

import {MMPA} from '../mmp-analysis/mmpa';
import {CLIFFS_TAB_TOOLTIP, FRAGMENTS_GRID_TOOLTIP, FRAGMENTS_TAB_TOOLTIP,
  GENERATIONS_TAB_TOOLTIP,
  MATCHED_MOLECULAR_PAIRS_TOOLTIP_CLIFFS, MATCHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS,
  MATCHED_MOLECULAR_PAIRS_TOOLTIP_TRANS, MMP_CONTEXT_PANE_CLASS, MMP_NAMES, SHOW_FRAGS_MODE,
  TrellisAxis, TrellisSortByProp, TrellisSortType} from './mmp-constants';

import {PaletteCodes, getPalette} from './palette';
import {getMmpTrellisPlot} from './mmp-frag-vs-frag';
import {getMmpScatterPlot, runMmpChemSpace, fillPairInfo} from './mmp-cliffs';
import {getGenerations} from './mmp-generations';

import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getEmbeddingColsNames} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {getMmpFilters, MmpFilters} from './mmp-filters';
import {getSigFigs} from '../../../utils/chem-common';
import {createLines} from './mmp-lines';
import {MmpPairedGrids} from './mmp-grids';
import {PackageFunctions} from '../../../package';
import {getZoomCoordinates} from '../../../utils/ui-utils';

export type MmpInput = {
  table: DG.DataFrame,
  molecules: DG.Column,
  activities: DG.ColumnList,
  fragmentCutoff: number
};

export type TrellisSorting = {
  [TrellisAxis.From]: SortType,
  [TrellisAxis.To]: SortType
}

export type SortType = {
  property: TrellisSortByProp,
  type?: TrellisSortType
}

export type SortData = {
  fragmentIdxs: number[],
  frequencies: number[],
  mw?: Float32Array,
}

type MMPHint = {
  text: string,
  element: HTMLElement,
  position: ui.hints.POSITION,
  class?: string,
  parentClass?: string,
}

export const STORAGE_NAME = 'mmpa';
export const KEY = 'hints';

export class MatchedMolecularPairsViewer extends DG.JsViewer {
  static TYPE: string = 'MMP';

  //properties
  molecules: string;
  moleculesColumnName: string | null = null;
  activities: string[] | null = null;
  diffTypes: string[] | null = null;
  scalings: string[] | null = null;
  fragmentCutoff: number | null;
  totalData: string;
  totalDataUpdated: boolean = false;
  onPropertyChangedObs : Subject<DG.Property | null> = new Subject<DG.Property | null>();
  hintsShown: {[key: string]: boolean} = {};

  moleculesCol: DG.Column | null = null;
  activitiesCols: DG.ColumnList | null = null;
  mmpa: MMPA | null = null;

  parentTable: DG.DataFrame | null = null;

  //mmpRules: MmpRules | null = null;
  colorPalette: PaletteCodes | null = null;

  pairedGrids: MmpPairedGrids | null = null;

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
  corrGrid: DG.Grid | null = null;
  generationsGridDiv = ui.div();

  rdkitModule: RDModule | null = null;
  currentTab = '';
  mmpFilters: MmpFilters | null = null;

  calculatedOnGPU: boolean | null = null;
  parentTableFilterBackup: DG.BitSet | null = null;
  mWCalulationsReady = false;

  fragSortingInfo: SortData = {fragmentIdxs: [], frequencies: []};
  sp: DG.ScatterPlotViewer | null = null;
  spAxesNames: string[] = [];
  tabs: DG.TabControl | null = null;
  mutationObserver: MutationObserver | null = null;

  tp: DG.Viewer | null = null;
  activityMeanNames: string[] = [];
  fragmentsDiv = ui.div();
  filterStatesUpdatedCondition: () => boolean = () => true;

  spCorrDiv = ui.divV([]);
  showFragmentsChoice: DG.InputBase | null = null;
  followCurrentRowInFragGrid: DG.InputBase | null = null;
  lastCurrentRowOnCliffsTab = - 1;
  lastOpenedHint: HTMLDivElement | null = null;
  showHints = true;

  constructor() {
    super();
    DG.debounce(this.onPropertyChangedObs, 1000).subscribe(this.onPropertyChangedDebounced.bind(this));
    //properties
    this.moleculesColumnName = this.addProperty('moleculesColumnName', DG.TYPE.COLUMN, '',
      {semType: DG.SEMTYPE.MOLECULE});
    this.activities = this.stringList('activities');
    this.fragmentCutoff = this.float('fragmentCutoff');
    this.molecules = this.string('molecules', '', {userEditable: false});

    this.totalData = this.string('totalData', 'null', {userEditable: false, includeInLayout: true});
    this.diffTypes = this.stringList('diffTypes', [], {userEditable: false, includeInLayout: true, nullable: false});
    this.scalings = this.stringList('scalings', [], {userEditable: false, includeInLayout: true, nullable: false});
  }

  onPropertyChangedDebounced() {
    if (!this.dataFrame)
      return;
    if (this.totalDataUpdated && this.moleculesColumnName) {
      this.moleculesCol = this.dataFrame.col(this.moleculesColumnName);
      this.activitiesCols = DG.DataFrame.fromColumns(this.dataFrame.columns.byNames(this.activities!)).columns;
      this.render();
      return;
    }

    this.moleculesCol = this.dataFrame.col(this.moleculesColumnName!);
    this.activitiesCols = DG.DataFrame.fromColumns(this.dataFrame.columns.byNames(this.activities!)).columns;
    if (this.moleculesColumnName && this.activities && this.fragmentCutoff && this.scalings) {
      this.render();
      return;
    }
  }

  onPropertyChanged(property: DG.Property | null): void {
    console.log(property!.name);
    super.onPropertyChanged(property);
    if (property?.name === 'totalData')
      this.totalDataUpdated = true;
    //for backward compatibility after changing molecules property to moleculesColumnName property
    if (property?.name === 'molecules' && property.get(this) !== '') {
      this.moleculesColumnName = property.get(this);
      this.molecules = '';
      return;
    }

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
        const errMsg = e instanceof Error ? e.message : e.toString();
        grok.log.error(e);
        grok.shell.error(errMsg);
      } finally {
        $(this.root).empty();
        if (this.tabs)
          this.root.appendChild(this.tabs.root);
        else
          this.close();
        progressMMP.close();
      }
    }
  }

  getTabs(): DG.TabControl {
    const decript1 = 'View all fragment substitutions found in the dataset';
    const decript2 = 'Analyze activity changes across fragment substitutions using a trellis plot';
    const decript3 = 'Molecule pairs analysis on 2d scatter plot';
    const decript4 = 'Generation of molecules based on obtained rules';

    const tabs = ui.tabControl(null, false);

    const transformationsTab = tabs.addPane(MMP_NAMES.TAB_TRANSFORMATIONS, () => {
      return this.getTransformationsTab();
    });
    ui.tooltip.bind(transformationsTab.header, decript1);

    const fragmentsTab = tabs.addPane(MMP_NAMES.TAB_FRAGMENTS, () => {
      this.prepareMwForSorting();
      //need timeout not to freeze when switching to tab
      setTimeout(() => {
        this.setupFragmentsTab();
      }, 100);
      ui.setUpdateIndicator(this.fragmentsDiv, true, 'Generating fragments trellis plot...');
      return this.fragmentsDiv;
    });
    ui.tooltip.bind(fragmentsTab.header, decript2);

    const cliffsTab = tabs.addPane(MMP_NAMES.TAB_CLIFFS, () => {
      return this.getCliffsTab();
    });
    ui.tooltip.bind(cliffsTab.header, decript3);

    const genTab = tabs.addPane(MMP_NAMES.TAB_GENERATION, () => {
      return this.getGenerationsTab();
    });
    ui.tooltip.bind(genTab.header, decript4);

    let firstCliffsOpening = true;
    tabs.onTabChanged.subscribe(() => {
      this.lastOpenedHint?.remove();
      const prevTab = this.currentTab;
      this.currentTab = tabs.currentPane.name;
      //when switching from cliffs tab - need to reset parent table filters
      if (prevTab === MMP_NAMES.TAB_CLIFFS && this.currentTab !== MMP_NAMES.TAB_CLIFFS)
        this.parentTable?.rows.requestFilter();

      //when switching from substitutions tab - hide pc plot
      if ((prevTab === MMP_NAMES.TAB_TRANSFORMATIONS || !prevTab) && this.currentTab !== MMP_NAMES.TAB_TRANSFORMATIONS)
        this.pairedGrids?.pcPlotContextPanelDiv.classList.add('mmpa-hide-context-panel');

      this.pairedGrids!.currentTab = tabs.currentPane.name as MMP_NAMES;

      if (tabs.currentPane.name == MMP_NAMES.TAB_TRANSFORMATIONS)
        this.pairedGrids!.switchToSubstitutionsTab();
      else if (tabs.currentPane.name == MMP_NAMES.TAB_FRAGMENTS) {
        this.pairedGrids!.switchToFragmentsTab();
        if (!fragmentsTab.content.classList.contains('mmpa-fragments-tab'))
          fragmentsTab.content.classList.add('mmpa-fragments-tab');
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_CLIFFS) {
        if (!cliffsTab.content.classList.contains('mmpa-cliffs-tab'))
          cliffsTab.content.classList.add('mmpa-cliffs-tab');

        this.sp!.root.append(this.mmpFilters!.filtersDiv);

        if (firstCliffsOpening)
          grok.shell.warning('Cutoff filters were applied for all activities');

        firstCliffsOpening = false;

        //requestFilter will call onRowsFiltering, where correct masks will be set
        this.parentTable?.rows.requestFilter();
        this.pairedGrids!.pairsGridCliffsTab.dataFrame.rows.requestFilter();
        //setting current arrow
        if (this.pairedGrids!.mmpGridTrans.dataFrame.currentRowIdx !== this.lastCurrentRowOnCliffsTab) {
          this.lastCurrentRowOnCliffsTab = this.pairedGrids!.mmpGridTrans.dataFrame.currentRowIdx;
          this.setCurrentArrow(this.lastCurrentRowOnCliffsTab);
        }

        //opening context panel with molecule pair info
        if (this.pairedGrids!.mmpGridTrans.dataFrame.currentRowIdx !== -1) {
          setTimeout(() => {
            grok.shell.windows.showContextPanel = true;
            grok.shell.o = fillPairInfo(this.mmpa!, this.pairedGrids!.mmpGridTrans.dataFrame.currentRowIdx,
              this.pairedGrids!.mmpGridTrans.dataFrame, this.parentTable!, this.rdkitModule!,
              this.moleculesCol!.name);
          }, 500);
        }
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_GENERATION) {
        grok.shell.windows.showContextPanel = true;
        grok.shell.o = this.spCorrDiv;
      }
    });

    return tabs;
  }

  getTransformationsTab(): HTMLElement {
    this.followCurrentRowInFragGrid = ui.input.bool('Follow current row', {
      value: false,
      onValueChanged: () => {
        this.pairedGrids!.followCurrentRowInFragmentsGrid = this.followCurrentRowInFragGrid!.value;
        this.pairedGrids!.refreshFragmentsAndPairs(false);
      },
    });

    ui.tooltip.bind(this.followCurrentRowInFragGrid.captionLabel,
      'If true, then \'Molecule pairs\' dataset below will be filtered according to current fragment pair');

    const mmPairsRoot1 = this.createGridDiv(MMP_NAMES.PAIRS_GRID,
      this.pairedGrids!.mmpGridTrans, MATCHED_MOLECULAR_PAIRS_TOOLTIP_TRANS,
      this.pairedGrids!.mmpGridTransMessage);

    mmPairsRoot1.prepend(
      ui.divText('No molecule pairs found. Try to change filter or select fragments pair from \'Fragments\' dataset.',
        'chem-mmpa-no-pairs-warning'));
    this.subs.push(this.pairedGrids!.showEmptyPairsWarningEvent.subscribe((showWarning: boolean) => {
      showWarning ? mmPairsRoot1.classList.add('chem-mmp-no-pairs') :
        mmPairsRoot1.classList.remove('chem-mmp-no-pairs');
    }));

    this.showFragmentsChoice = ui.input.choice('', {items: [SHOW_FRAGS_MODE.All, SHOW_FRAGS_MODE.Current],
      nullable: false, value: SHOW_FRAGS_MODE.All,
      onValueChanged: (value) => {
        this.pairedGrids!.fragsShowAllMode = value === SHOW_FRAGS_MODE.All;
        if (value === SHOW_FRAGS_MODE.All) {
          this.pairedGrids!.fpMaskByMolecule.setAll(true);
          this.pairedGrids!.mmpMaskTrans.setAll(true);
          this.pairedGrids!.fpGrid.dataFrame.rows.requestFilter();
          this.pairedGrids!.mmpGridTrans.dataFrame.rows.requestFilter();
        } else
          this.pairedGrids!.refreshFragmentsAndPairs(true);
      }});
    this.showFragmentsChoice.root.classList.add('chem-mmp-fragments-grid-mode-choice');
    ui.tooltip.bind(this.showFragmentsChoice.input,
      `Select whether to show all fragments pairs or only pairs for current molecule in the initial dataset`);


    this.followCurrentRowInFragGrid.classList.add('chem-mmp-fragments-grid-follow-current-row');

    const fpGrid = this.createGridDiv(MMP_NAMES.FRAGMENTS_GRID,
      this.pairedGrids!.fpGrid, FRAGMENTS_GRID_TOOLTIP, this.pairedGrids!.fpGridMessage,
      ui.divH([this.showFragmentsChoice.root, this.followCurrentRowInFragGrid.root]));
    fpGrid.prepend(
      // eslint-disable-next-line max-len
      ui.divText('No substitutions found. Try to change filters or select another molecule if you are in \'Current\' mode.',
        'chem-mmpa-no-fragments-error'));
    this.subs.push(this.pairedGrids!.showErrorEvent.subscribe((showError: boolean) => {
      showError ? fpGrid.classList.add('chem-mmp-no-fragments') : fpGrid.classList.remove('chem-mmp-no-fragments');
    }));

    const hints: MMPHint[] = [
      {
        element: fpGrid,
        text: FRAGMENTS_GRID_TOOLTIP,
        position: ui.hints.POSITION.LEFT,
        class: 'chem-mmp-active-hint-element-horz',
      },
      {
        element: mmPairsRoot1,
        text: MATCHED_MOLECULAR_PAIRS_TOOLTIP_TRANS,
        position: ui.hints.POSITION.LEFT,
        class: 'chem-mmp-active-hint-element-horz',
      },
      {
        element: grok.shell.tv.grid.root,
        text: `Observe two molecules from the selected pair pinned on top of the grid.`,
        position: ui.hints.POSITION.RIGHT,
        class: 'chem-mmp-active-hint-element-horz',
      },
      {
        element: this.showFragmentsChoice.root,
        // eslint-disable-next-line max-len
        text: `Change mode to 'Current molecule' to filter all subtitutions for current molecule from the initial dataset on the left.`,
        position: ui.hints.POSITION.LEFT,
        parentClass: 'chem-mmp-active-hint-adjust-vert-2px',
      },
      {
        element: grok.shell.tv.grid.root,
        text: `Change current row to see the changes in Fragments grid`,
        position: ui.hints.POSITION.RIGHT,
        class: 'chem-mmp-active-hint-element-horz',
      },
      {
        element: this.followCurrentRowInFragGrid.root,
        // eslint-disable-next-line max-len
        text: `Change to true and click any row in 'Fragments' grid to filter molecule pairs with corresponding substitution in 'Molecule pairs' grid`,
        position: ui.hints.POSITION.RIGHT,
        class: 'chem-mmp-active-hint-element-horz',
      },
    ];

    setTimeout(() => {
      this.setupHint(hints, 0);
    }, 1000);

    const gridsDiv = ui.splitV([
      fpGrid,
      mmPairsRoot1,
    ], {}, true);

    return gridsDiv;
  }

  setupHint(hints: MMPHint[], i: number) {
    if (!this.showHints)
      return;
    const doNotShow = ui.link('Do not show', () => {
      this.lastOpenedHint?.remove();
      this.showHints = false;
      grok.userSettings.add(STORAGE_NAME, KEY, 'false');
    });
    doNotShow.style.paddingTop = '11px';
    const addHint = (i: number, el: HTMLElement) => {
      if (!this.showHints)
        return;
      this.lastOpenedHint = ui.hints.addHint(hints[i].element, el, hints[i].position);
      this.lastOpenedHint.classList.add(`chem-mmp-hint-${i}`);
      hints[i].element.classList.add('chem-mmp-active-hint-element');
      if (hints[i].class)
        hints[i].element.classList.add(hints[i].class);
      if (hints[i].parentClass)
        hints[i].element.parentElement?.classList.add(hints[i].parentClass);
    };
    this.lastOpenedHint?.remove();
    const hintContent = ui.div();
    hintContent.append(ui.divText(hints[i].text));
    const nextButton = i !== hints.length - 1 ? ui.button('Next', () => this.setupHint(hints, i + 1)) :
      ui.button('Close', () => this.lastOpenedHint?.remove());
    const buttons = ui.divH([], {style: {justifyContent: 'space-between'}});
    buttons.append(nextButton);
    buttons.append(doNotShow);
    hintContent.append(buttons);
    addHint(i, hintContent);
    const hintMutationObserver = new MutationObserver((mutationsList) => {
      for (let j = 0; j < mutationsList.length; j++) {
        for (const node of Array.from(mutationsList[j].removedNodes)) {
          const removedHint = (node as HTMLElement).classList.contains(`chem-mmp-hint-${i}`);
          if (removedHint) {
            hints[i].element.classList.remove('chem-mmp-active-hint-element');
            if (hints[i].class)
              hints[i].element.classList.remove(hints[i].class);
            if (hints[i].parentClass)
              hints[i].element.parentElement?.classList.remove(hints[i].parentClass);
            hintMutationObserver?.disconnect();
          }
        }
      }
    });
    hintMutationObserver.observe(document.body, {attributes: true, childList: true});
  }

  setupFragmentsTab(): void {
    //Fragments tab
    this.tp = getMmpTrellisPlot(this.pairedGrids!.fpGrid, this.activityMeanNames, this.colorPalette!);
    this.tp.onEvent('d4-trellis-plot-inner-viewer-clicked').subscribe((cats) => {
      //grok.shell.info(cats);
      this.pairedGrids?.refilterMatchedPairsByFragments(cats);
      this.pairedGrids?.mmpGridFrag.dataFrame.rows.requestFilter();
    });

    //setup trellis filters
    const trellisTv = DG.TableView.create(this.pairedGrids!.fpGrid.dataFrame, false);
    this.pairedGrids!.filters = trellisTv.getFiltersGroup();

    const trellisHeader = ui.h1('Fragment vs Fragment', 'chem-mmpa-transformation-tab-header');

    const filterIcon = ui.icons.filter(() => {
      ui.showPopup(ui.div(this.pairedGrids!.filters!.root, 'chem-mmp-trellis-filters-div'),
        filterIcon, {vertical: true});
    }, 'Fragments filters');
    filterIcon.classList.add('chem-mmpa-fragments-filters-icon');

    const trellisSortState: TrellisSorting = {
      [TrellisAxis.From]: {property: TrellisSortByProp.Frequency, type: TrellisSortType.Desc},
      [TrellisAxis.To]: {property: TrellisSortByProp.Frequency, type: TrellisSortType.Desc},
    };

    const summaryColsButton = ui.div('', 'mmp-trellis-summary-column');
    ui.tooltip.bind(summaryColsButton, 'Select columns to show in trellis plot');

    //create trellis legend
    const trellisLegend = ui.divH([], 'mmpa-trellis-legend');
    this.updateTrellisLegend(trellisLegend, this.tp!.getOptions().look.innerViewerLook.columnNames);

    this.tp.root.prepend(trellisHeader);
    const tpDiv = ui.splitV([
      ui.box(
        ui.divH([
          trellisHeader, filterIcon, summaryColsButton,
          this.helpButton('chem-mmpa-grid-help-icon', FRAGMENTS_TAB_TOOLTIP),
          trellisLegend,
        ], {style: {width: '100%'}}),
        {style: {maxHeight: '30px'}},
      ),
      this.tp.root,
    ], {style: {width: '100%', height: '100%'}});

    this.tp.onEvent('d4-viewer-rendered').subscribe(() => {
      this.updateTrellisLegend(trellisLegend, this.tp!.getOptions().look.innerViewerLook.columnNames);
      this.createSortIcon(trellisSortState, TrellisAxis.From, this.tp!, 'chem-mmpa-fragments-sort-icon-x-axis');
      this.createSortIcon(trellisSortState, TrellisAxis.To, this.tp!, 'chem-mmpa-fragments-sort-icon-y-axis');
      const tpButtons = Array.from(this.tp!.root.getElementsByTagName('button'));
      if (tpButtons.length) {
        //workaround to disable aggregation functions selections
        tpButtons[0].onmousedown = () => {
          if (!this.mutationObserver) {
            this.mutationObserver = new MutationObserver((mutationsList) => {
              for (let i = 0; i < mutationsList.length; i++) {
                for (const node of Array.from(mutationsList[i].addedNodes)) {
                  const dlgHeaders = Array.from((node as HTMLElement).getElementsByClassName('d4-dialog-title'))
                    .filter((el) => (el as HTMLElement).innerText === 'Edit columns aggregations');
                  if (dlgHeaders.length) {
                    dlgHeaders[0].parentElement?.parentElement?.classList.add('mmp-trellis-summary-col-dlg');
                    this.mutationObserver?.disconnect();
                    this.mutationObserver = null;
                    return;
                  }
                }
              }
            });
            this.mutationObserver.observe(document.body, {attributes: true, childList: true});
          }
        };
        ui.empty(summaryColsButton);
        summaryColsButton.append(tpButtons[0]);
      }
    });

    this.sortTrellis(TrellisAxis.From, trellisSortState[TrellisAxis.From], this.tp);
    this.sortTrellis(TrellisAxis.To, trellisSortState[TrellisAxis.To], this.tp);

    const mmPairsRoot2 = this.createGridDiv(MMP_NAMES.PAIRS_GRID,
      this.pairedGrids!.mmpGridFrag, MATCHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS, this.pairedGrids!.mmpGridFragMessage);

    const hints: MMPHint[] = [
      {
        element: this.tp.root,
        text: FRAGMENTS_TAB_TOOLTIP,
        position: ui.hints.POSITION.LEFT,
      },
      {
        element: filterIcon,
        text: `Use trellis plot filters to filter fragments of interest`,
        position: ui.hints.POSITION.LEFT,
      },
    ];

    setTimeout(() => {
      this.setupHint(hints, 0);
    }, 1000);

    ui.empty(this.fragmentsDiv);
    ui.setUpdateIndicator(this.fragmentsDiv, false);
    this.fragmentsDiv.append(ui.splitV([
      tpDiv,
      mmPairsRoot2,
    ], {style: {width: '100%', height: '100%'}}, true));
  }

  updateTrellisLegend(trellisLegend: HTMLDivElement, selectedActivities: string[]) {
    ui.empty(trellisLegend);
    for (let i = 0; i < selectedActivities.length; i++) {
      const div = ui.divText(selectedActivities[i], {style: {color: this.colorPalette!.hex[i]}});
      div.classList.add('mmpa-trellis-legend-item');
      //show tooltip only in case name is truncated
      ui.tooltip.bind(div, () => div.scrollWidth > div.clientWidth ? selectedActivities[i] : null);
      trellisLegend.append(div);
    }
  }

  getCliffsTab(): HTMLElement {
    const {linesIdxs, lines, linesActivityCorrespondance} = createLines(this.mmpa!, this.colorPalette!);
    this.linesIdxs = linesIdxs;
    this.lines = lines;
    this.spAxesNames = getEmbeddingColsNames(this.parentTable!).map((it) => `~${it}`);
    this.linesActivityCorrespondance = linesActivityCorrespondance;
    this.linesMask = new BitArray(linesIdxs.length);

    this.mmpFilters = getMmpFilters(this.activities!, this.mmpa!.allCasesBased.maxActs);
    this.cutoffMasks = new Array<DG.BitSet>(this.mmpFilters.activitySliderInputs.length);
    this.totalCutoffMask = DG.BitSet.create(this.parentTable!.rowCount);
    this.setupFilters(this.mmpFilters, linesActivityCorrespondance);
    console.log(`created mmpa filters`);

    this.sp = getMmpScatterPlot(this.parentTable!, this.spAxesNames,
      this.moleculesCol!.name, this.activitiesCols!.byIndex(0).name);
    //show scatter plot context menu instead of mmp viewer's context menu
    this.subs.push(grok.events.onContextMenu.subscribe((e) => {
      if (e.causedBy && e.causedBy.target && this.sp!.root.contains(e.causedBy.target)) {
        e.causedBy.preventDefault();
        e.causedBy.stopPropagation();
        e.causedBy.stopImmediatePropagation();
      }
    }));

    const [linesEditor, chemSpaceParams] = runMmpChemSpace(this.parentTable!, this.moleculesCol!, this.sp, lines,
      linesIdxs, linesActivityCorrespondance, this.pairedGrids!.mmpGridTrans.dataFrame, this.mmpa!, this.rdkitModule!,
      this.spAxesNames);

    const progressBarSpace = DG.TaskBarProgressIndicator.create(`Running Chemical space...`);
    this.mmpa!.chemSpace(chemSpaceParams).then((res) => {
      const embeddings = res.coordinates;
      for (const col of embeddings)
        this.parentTable!.columns.replace(col.name, col);
      //workaround for case when sp is opened for the first time with minimal height
      spDiv.style.height = '800px';
      this.totalData = this.mmpa!.toJSON();
      progressBarSpace.close();
    }).catch((error: any) => {
      const errorStr = `Scatter plot haven't been completed due to error: ${error}`;
      grok.shell.error(errorStr);
    }).finally(() => ui.setUpdateIndicator(this.sp!.root, false));
    this.sp.root.style.width = '100%';

    this.totalCutoffMask!.setAll(true);
    this.linesMask!.setAll(true);

    this.linesRenderer = linesEditor;

    this.linesRenderer!.lineClicked.subscribe((event: MouseOverLineEvent) => {
      this.linesRenderer!.currentLineId = event.id;
      if (event.id !== -1) {
        const pairId = this.linesIdxs![event.id];
        this.pairedGrids!.pairsGridCliffsTab.dataFrame.currentRowIdx = pairId;
      }
    });

    this.subs.push(this.parentTable!.onRowsFiltering.subscribe(() => {
      if (this.currentTab === MMP_NAMES.TAB_CLIFFS)
        this.parentTable!.filter.and(this.totalCutoffMask!);
    }));

    this.refilterCliffs(this.mmpFilters.activitySliderInputs.map((si) => si.value),
      this.mmpFilters.activityActiveInputs.map((ai) => ai.value));

    const mmPairsRoot3 = this.createGridDiv(MMP_NAMES.PAIRS_GRID, this.pairedGrids!.pairsGridCliffsTab,
      MATCHED_MOLECULAR_PAIRS_TOOLTIP_CLIFFS, this.pairedGrids!.pairsGridCliffsTabMessage);
    mmPairsRoot3.classList.add('mmp-pairs-grid-cliffs-tab', 'cliffs-opened');

    this.pairedGrids!.pairsGridCliffsTab.dataFrame.onCurrentRowChanged.subscribe(() => {
      const currentRowIdx = this.pairedGrids!.pairsGridCliffsTab.dataFrame.currentRowIdx;
      if (currentRowIdx !== -1 && this.currentTab === MMP_NAMES.TAB_CLIFFS) {
        this.lastCurrentRowOnCliffsTab = currentRowIdx;
        this.setCurrentArrow(this.lastCurrentRowOnCliffsTab);
      }
    });

    let pairsOpened = true;
    const cliffsHeader = ui.h1('Chemical space', 'chem-mmpa-transformation-tab-header');
    const showPairs = ui.input.bool('Show pairs', {value: pairsOpened, onValueChanged: () => {
      pairsOpened = !pairsOpened;
      if (pairsOpened)
        mmPairsRoot3.classList.replace('cliffs-closed', 'cliffs-opened');
      else
        mmPairsRoot3.classList.replace('cliffs-opened', 'cliffs-closed');
    }});
    showPairs.root.classList.add('chem-mmp-show-pairs-checkbox');

    ui.setUpdateIndicator(this.sp.root, true, 'Genarating cliffs scatter plot...');

    const hints: MMPHint[] = [
      {
        element: this.sp.root,
        // eslint-disable-next-line max-len
        text: `Scatter plot where similar molecules are close to each other. Lines connect matched molecular pairs. Arrow points to a molecule with a greater activity value. Click on a line to show molecule pair in the table below and show details in a context panel.`,
        position: ui.hints.POSITION.LEFT,
        class: 'chem-mmp-active-hint-element-horz',
      },
      {
        element: this.mmpFilters.filtersDiv,
        // eslint-disable-next-line max-len
        text: `Use filters on the scatter plot to change activity difference cutoff or switch of/on the lines for any activity.`,
        position: ui.hints.POSITION.LEFT,
      },
      {
        element: mmPairsRoot3,
        text: MATCHED_MOLECULAR_PAIRS_TOOLTIP_CLIFFS,
        position: ui.hints.POSITION.LEFT,
        class: 'chem-mmp-active-hint-element-horz',
      },
      {
        element: showPairs.root,
        text: `Use checkbox to show/hide molecule pairs grid`,
        position: ui.hints.POSITION.LEFT,
        parentClass: 'chem-mmp-active-hint-adjust-vert-5px',
      },
    ];

    setTimeout(() => {
      this.setupHint(hints, 0);
    }, 1000);

    const spDiv = ui.splitV([
      ui.box(
        ui.divH([cliffsHeader, showPairs.root, this.helpButton('chem-mmpa-grid-help-icon', CLIFFS_TAB_TOOLTIP)]),
        {style: {maxHeight: '30px'}},
      ),
      this.sp.root,
    ]);

    return ui.splitV([
      spDiv,
      mmPairsRoot3,
    ], {style: {width: '100%', height: '100%'}}, true);
  }

  setCurrentArrow(currentRowIdx: number) {
    const fromIdx = this.pairedGrids!.pairsGridCliffsTab.dataFrame.get(MMP_NAMES.PAIRNUM_FROM, currentRowIdx);
    const toIdx = this.pairedGrids!.pairsGridCliffsTab.dataFrame.get(MMP_NAMES.PAIRNUM_TO, currentRowIdx);
    //this.lineIdxs contain idxs if pairs from pairs dataset
    let currentLineIdx: number | null = null;
    for (let i = 0; i < this.linesIdxs!.length; i++) {
      if (this.linesIdxs![i] === currentRowIdx && this.linesRenderer!.visibility?.getBit(i)) {
        currentLineIdx = i;
        break;
      }
    }
    if (currentLineIdx) {
      this.linesRenderer!.currentLineId = currentLineIdx;
      const {zoomLeft, zoomRight, zoomTop, zoomBottom} = getZoomCoordinates(
          this.sp!.viewport.width,
          this.sp!.viewport.height,
          this.sp!.dataFrame.get(this.spAxesNames[0], fromIdx),
          this.sp!.dataFrame.get(this.spAxesNames[1], fromIdx),
          this.sp!.dataFrame.get(this.spAxesNames[0], toIdx),
          this.sp!.dataFrame.get(this.spAxesNames[1], toIdx),
      );
        this.sp!.zoom(zoomLeft,
          zoomTop,
          zoomRight,
          zoomBottom);
    }
  }

  //cliffs filters on scatter plot
  setupFilters(mmpFilters: MmpFilters, linesActivityCorrespondance: Uint32Array): void {
    for (let i = 0; i < mmpFilters.activitySliderInputs.length; i++) {
      DG.debounce(mmpFilters.activityActiveInputs[i].onChanged, 300).subscribe(() => {
        this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
          mmpFilters.activityActiveInputs.map((ai) => ai.value));
        this.parentTable?.rows.requestFilter();
      });

      DG.debounce(mmpFilters.activitySliderInputs[i].onChanged, 300).subscribe(() => {
        mmpFilters.activityValuesDivs[i].innerText = mmpFilters.activitySliderInputs[i].value === 0 ? '0' :
          getSigFigs(mmpFilters.activitySliderInputs[i].value, 4).toString();
        this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
          mmpFilters.activityActiveInputs.map((ai) => ai.value));
        this.parentTable?.rows.requestFilter();
      });

      mmpFilters.activityColorInputs[i].value = this.colorPalette!.hex[i];
      mmpFilters.activityColorInputs[i].onChanged.subscribe(() => {
        const progressRendering = DG.TaskBarProgressIndicator.create(`Changing colors...`);

        //refresh lines
        for (let i =0; i < mmpFilters.activityColorInputs.length; i++) {
          this.colorPalette!.hex[i] = mmpFilters.activityColorInputs[i].value;
          this.colorPalette!.numerical[i] = DG.Color.fromHtml(mmpFilters.activityColorInputs[i].value);
          this.colorPalette!.rgb[i] = DG.Color.toRgb(this.colorPalette!.numerical[i]);
          this.colorPalette!.rgbCut[i] = this.colorPalette!.rgb[i].replace('rgb(', '').replace(')', '');
        }

        const colors = this.lines!.colors;
        for (let i = 0; i < colors!.length; i++)
          colors![i] = this.colorPalette!.rgbCut[linesActivityCorrespondance[i]];

        const lines: ILineSeries = {
          from: this.lines!.from,
          to: this.lines!.to,
          drawArrows: true,
          colors: colors,
          arrowSize: 10,
          skipMultiLineCalculation: true,
          width: 0.5,
        };

        this.lines = lines;
        this.linesRenderer!.linesToRender = lines;

        //refresh trellis plot
        const schemes = new Array<any>(mmpFilters.activityColorInputs.length);
        for (let i = 0; i < mmpFilters.activityColorInputs.length; i++)
          schemes[i] = [this.colorPalette!.numerical[i]];

        this.tp?.setOptions({'innerViewerLook': {'colorSchemes': schemes}});
        progressRendering.close();
      });

      this.cutoffMasks![i] = DG.BitSet.create(this.parentTable!.rowCount);
      this.cutoffMasks![i].setAll(true);
    }
  }

  getGenerationsTab(): HTMLElement {
    getGenerations(this.mmpa!, this.pairedGrids!.fpGrid).then(([genGrid, corrGrid]) => {
      this.generationsGrid = genGrid;
      this.corrGrid = corrGrid;

      this.corrGrid.dataFrame.meta.formulaLines.addLine({
        title: 'Identity',
        formula: '${Observed} = ${Predicted}',
        color: '#717581',
        width: 1,
        visible: true,
      });

      ui.setUpdateIndicator(this.generationsGridDiv, false);
      ui.empty(this.generationsGridDiv);
      this.generationsGridDiv.append(
        this.createGridDiv('Generated Molecules', this.generationsGrid!, GENERATIONS_TAB_TOOLTIP, ui.div()));

      for (const activity of this.activities!) {
        const generationsSp = DG.Viewer.scatterPlot(this.corrGrid?.dataFrame!, {
          x: 'Observed',
          y: 'Predicted',
          filter: `\${${MMP_NAMES.ACTIVITY}} == "${activity}" `,
          zoomAndFilter: 'no action',
          color: 'Activity',
          showXSelector: true,
          showXAxis: true,
          showYSelector: true,
          showYAxis: true,
          showColorSelector: true,
          showSizeSelector: true,
          markerDefaultSize: this.mmpa!.initData.molecules.length > 10000 ? 1 : 2,
          markerType: 'circle',
          showRegressionLine: true,
          legendVisibility: 'Never',
        });
        const header = ui.h2(`${activity}: Observed vs Predicted`, 'chem-mmpa-generation-tab-cp-header');
        generationsSp.root.prepend(header);
        const spActivityCorrDiv = ui.div(ui.splitV([
          ui.box(
            ui.divH([header]),
            {style: {maxHeight: '30px'}},
          ),
          generationsSp.root,
        ]));
        this.spCorrDiv.append(spActivityCorrDiv);
      }
      ui.setUpdateIndicator(this.spCorrDiv, false);
      this.spCorrDiv.classList.add(MMP_CONTEXT_PANE_CLASS);
      grok.shell.windows.showContextPanel = true;
      grok.shell.o = this.spCorrDiv;

      const hints: MMPHint[] = [
        {
          element: this.generationsGridDiv,
          text: GENERATIONS_TAB_TOOLTIP,
          position: ui.hints.POSITION.LEFT,
          class: 'chem-mmp-active-hint-element-horz',
        },
      ];

      setTimeout(() => {
        this.setupHint(hints, 0);
      }, 1000);
    }).catch((error: any) => {
      const errorStr = `Generations haven't been completed due to error: ${error}`;
      ui.setUpdateIndicator(this.generationsGridDiv, false);
      ui.setUpdateIndicator(this.spCorrDiv, false);
      ui.empty(this.generationsGridDiv);
      this.generationsGridDiv.append(ui.divText(errorStr));
      grok.shell.error(errorStr);
    });
    ui.setUpdateIndicator(this.generationsGridDiv, true, 'Generations in progress...');
    ui.setUpdateIndicator(this.spCorrDiv, true, 'Generating correlation scatter plot...');
    return this.generationsGridDiv;
  }

  createGridDiv(name: string, grid: DG.Grid, helpTooltip: string, messageBox: HTMLElement, extraEl?: HTMLElement) {
    const header = ui.h1(name, 'chem-mmpa-transformation-tab-header');
    const helpBUtton = this.helpButton('chem-mmpa-grid-help-icon', helpTooltip);
    grid.root.prepend(header);
    return ui.splitV([
      ui.box(
        ui.divH([ui.divH([header, extraEl ?? ui.div(),
          this.addToWorkspaceButton(grid.dataFrame, name, 'chem-mmpa-add-to-workspace-button'), helpBUtton]),
        messageBox], {style: {justifyContent: 'space-between'}}),
        {style: {maxHeight: '32px'}},
      ),
      grid.root,
    ], {style: {width: '100%', height: '100%'}});
  };

  addToWorkspaceButton(table: DG.DataFrame, name: string, className: string) {
    const button = ui.icons.add(() => {
      const clonedTable = table.clone();
      clonedTable.name = name;
      const tv = grok.shell.addTableView(clonedTable);
      //without setTimeout tableView is not set as current view
      setTimeout(() => grok.shell.v = tv, 10);
    }, 'Add table to workspace');
    button.classList.add(className);
    return button;
  };

  helpButton(className: string, tooltip: string) {
    const button = ui.icons.help(() => {
      grok.shell.windows.showHelp = true;
      grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#matched-molecular-pairs');
    });
    ui.tooltip.bind(button, () => ui.divText(tooltip, {style: {width: '300px'}}));
    button.classList.add(className);
    return button;
  };

  createSortIcon(trellisSortState: TrellisSorting, axis: TrellisAxis, tp: DG.Viewer, className: string) {
    const tpAxisSelectorArr = Array.from(tp.root.getElementsByClassName('d4-column-selector-column'))
      .filter((it) => (it as HTMLElement).innerText === axis);
    if (tpAxisSelectorArr.length) {
      const tpAxisSelector = tpAxisSelectorArr[0] as HTMLElement;
      const sortIconExists = tpAxisSelector.parentElement!.parentElement!.getElementsByClassName(className).length;
      if (!sortIconExists) {
        const sortPropChoice = ui.input.choice('Sort by', {
          value: trellisSortState[axis].property, items: Object.keys(TrellisSortByProp), nullable: false,
          onValueChanged: () => {
            if (trellisSortState[axis].property !== sortPropChoice.value) {
              if (sortPropChoice.value === TrellisSortByProp.None)
                sortTypeChoice.root.classList.add('chem-mmp-trellis-sort-type-invisible');
              else
                sortTypeChoice.root.classList.remove('chem-mmp-trellis-sort-type-invisible');
              trellisSortState[axis].property = sortPropChoice.value as TrellisSortByProp;
              this.sortTrellis(axis, trellisSortState[axis], tp);
            }
          },
        });
        const sortTypeChoice = ui.input.radio('', {
          value: trellisSortState[axis].type, items: Object.keys(TrellisSortType), nullable: false,
          onValueChanged: () => {
            if (trellisSortState[axis].type !== sortTypeChoice.value) {
              trellisSortState[axis].type = sortTypeChoice.value as TrellisSortType;
              this.sortTrellis(axis, trellisSortState[axis], tp);
            }
          },
        });
        if (trellisSortState[axis].property === TrellisSortByProp.None)
          sortTypeChoice.root.classList.add('chem-mmp-trellis-sort-type-invisible');
        const sortIcon = ui.iconFA('sort-alt', () => {
          ui.showPopup(ui.inputs([
            sortPropChoice, sortTypeChoice,
          ], 'chem-mmp-trellis-plot-sort-div'), sortIcon);
        }, `Sort ${axis} axis`);
        sortIcon.classList.add(className);
          tpAxisSelector.parentElement!.parentElement!.append(sortIcon);
          tpAxisSelector.onclick = (e: MouseEvent) => {
            e.stopImmediatePropagation();
            e.preventDefault();
          };
      }
    }
  }

  sortTrellis(axis: TrellisAxis, sorting: SortType, tp: DG.Viewer) {
    //sorting will reset the filter so need to save it for backup
    const filterBackup = this.pairedGrids!.fpGrid.dataFrame.filter.clone();
    switch (sorting.property) {
    case TrellisSortByProp.None:
      this.sortTrellisByProp(axis, tp);
      break;
    case TrellisSortByProp.Frequency:
      this.sortTrellisByProp(axis, tp, this.fragSortingInfo.frequencies,
        sorting.type);
      break;
    case TrellisSortByProp.MW:
      if (!this.mWCalulationsReady) {
        grok.shell.warning('MW calculations for fragments in progress. Please try again later');
        return;
      } else if (!this.fragSortingInfo.mw || this.fragSortingInfo.mw.every((it) => it === 0)) {
        grok.shell.warning('Molecular weights haven\'t been calculated');
        return;
      }
      this.sortTrellisByProp(axis, tp, this.fragSortingInfo.mw, sorting.type);
    }
    this.pairedGrids!.fpGrid.dataFrame.filter.copyFrom(filterBackup);
  }

  sortTrellisByProp(axis: TrellisAxis, tp: DG.Viewer, sortData?: number[] | Float32Array,
    type?: TrellisSortType) {
    const fragCol = this.pairedGrids!.fpGrid.dataFrame.col(axis);
    const axisName = tp.props.yColumnNames[0] === axis ? 'yColumnNames' :
      tp.props.xColumnNames[0] === axis ? 'xColumnNames' : null;

    const frags = this.fragSortingInfo.fragmentIdxs.map((idx) => this.mmpa?.frags.idToName[idx]) as string[];
    const indexesArray = [...Array(this.fragSortingInfo.fragmentIdxs.length).keys()];

    if (fragCol) {
      let cats: string[] = [];
      if (!type)
        cats = fragCol.categories.map((it) => it).sort();
      else {
        if (this.fragSortingInfo.fragmentIdxs && sortData) {
          switch (type) {
          case TrellisSortType.Asc:
            cats = indexesArray.sort((a: number, b: number) => sortData[a] - sortData[b] > 0 ? 1 : -1)
              .map((it) => frags[it]);
            break;
          case TrellisSortType.Desc:
            cats = indexesArray.sort((a: number, b: number) => sortData[a] - sortData[b] > 0 ? -1 : 1)
              .map((it) => frags[it]);
            break;
          }
        }
      }
      if (cats.length)
        fragCol.setCategoryOrder(cats);
    }
    if (axisName)
      tp.props[axisName] = tp.props[axisName];
  }

  fillAll(mmpInput: MmpInput, palette: PaletteCodes, mmpa: MMPA, diffs: Array<Float32Array>,
    pairedGrids: MmpPairedGrids, activityMeanNames: string[]): void {
    this.rdkitModule = getRdKitModule();
    this.activityMeanNames = activityMeanNames;
    this.parentTable = mmpInput.table;
    this.diffs = diffs;
    this.colorPalette = palette;
    this.mmpa = mmpa;


    //main grids
    this.pairedGrids = pairedGrids;
    this.pairedGrids.setupGrids();

    //tabs
    this.tabs = this.getTabs();

    this.root.append(this.tabs.root);

    this.calculatedOnGPU = mmpa.gpu;
  }

  async runMMP(mmpInput: MmpInput) {
    //console.profile('MMP');
    const showHints = grok.userSettings.getValue(STORAGE_NAME, KEY);
    if (showHints === 'false')
      this.showHints = false;
    const moleculesArray = mmpInput.molecules.toList();
    const variates = mmpInput.activities.length;
    const activitiesArrays = new Array<Float32Array>(variates);
    const activitiesNames = new Array<string>(variates);
    for (let i = 0; i < variates; i++) {
      activitiesArrays[i] = mmpInput.activities.byIndex(i).asDoubleList();
      activitiesNames[i] = mmpInput.activities.byIndex(i).name;
    }

    let mmpa: MMPA;

    try {
      if (this.totalDataUpdated) {
        mmpa = await MMPA.fromData(
          mmpInput.molecules.name, this.totalData, moleculesArray,
          activitiesArrays, activitiesNames, this.diffTypes!, this.fragSortingInfo);
        this.totalDataUpdated = false;
      } else {
        mmpa = await MMPA.init(
          mmpInput.molecules.name, moleculesArray, mmpInput.fragmentCutoff,
          activitiesArrays, activitiesNames, this.diffTypes!, this.fragSortingInfo);
      }
    } catch (err: any) {
      const errMsg = err instanceof Error ? err.message : err.toString();
      grok.log.error(err);
      grok.shell.error(errMsg);
      throw new Error(errMsg);
    }

    const palette = getPalette(mmpInput.activities.length);
    const activityMeanNames = Array<string>(mmpa.initData.activitiesCount);
    for (let i = 0; i < mmpa.initData.activitiesCount; i++) {
      const name = MMP_NAMES.MEANDIFF + ' ' + mmpa.initData.activitiesNames[i];
      activityMeanNames[i] = name;
    }

    //Transformations tab
    const pairedGrids = new MmpPairedGrids(this.subs, mmpInput, mmpa, activityMeanNames, palette);


    this.fillAll(mmpInput, palette, mmpa, mmpa.allCasesBased.diffs, pairedGrids, activityMeanNames);

    this.totalData = mmpa.toJSON();
    //console.profileEnd('MMP');
  }


  refilterCliffs(cutoffs: number[], isActiveVar: boolean[]): void {
    for (let i = 0; i < this.cutoffMasks!.length; i++)
      this.cutoffMasks![i].setAll(false);

    this.totalCutoffMask!.setAll(false);
    this.linesMask!.setAll(false);
    this.pairedGrids?.pairsMaskCliffsTab.setAll(false);

    //setting activity associated masks
    for (let i = 0; i < this.lines!.from.length; i++) {
      const activityNumber = this.linesActivityCorrespondance![i];
      //line is idx of pair in pairs dataset
      const line = this.linesIdxs![i];
      //if 'isActive' checkbox for variable is unset
      //setting line to invisible and continue, otherwise look at cutoff
      if (!isActiveVar[activityNumber]) {
        this.linesMask!.setBit(i, false, false);
        continue;
      }
      if (this.diffs![activityNumber][line] >= cutoffs[activityNumber]) {
        this.cutoffMasks![activityNumber].set(this.lines!.from[i], true, false);
        this.cutoffMasks![activityNumber].set(this.lines!.to[i], true, false);
        this.linesMask!.setBit(i, true, false);
        this.pairedGrids?.pairsMaskCliffsTab.set(line, true);
      }
    }

    //setting resulting masks
    //TODO: refine
    for (let j = 0; j < this.totalCutoffMask!.length; j++) {
      let res = false;
      for (let i = 0; i < this.cutoffMasks!.length; i++)
        res = res || this.cutoffMasks![i].get(j);
      this.totalCutoffMask!.set(j, res, false);
    }
    this.linesRenderer!.linesVisibility = this.linesMask!;
  }

  async prepareMwForSorting() {
    const frags = this.fragSortingInfo.fragmentIdxs.map((idx) => this.mmpa?.frags.idToName[idx]) as string[];
    PackageFunctions.getDescriptors(DG.Column.fromStrings('smiles', frags), ['MolWt']).then((res: DG.DataFrame) => {
      this.fragSortingInfo.mw = new Float32Array(this.fragSortingInfo.fragmentIdxs.length);
      let errorCount = 0;
      frags.forEach((key: string, idx) => {
        let resMW = res.get('MolWt', idx);
        if (!resMW || typeof resMW === 'string') {
          errorCount++;
          resMW = undefined;
        }
        this.fragSortingInfo.mw![idx] = resMW ?? 0;
      });
      this.mWCalulationsReady = true;
      if (errorCount > 0)
        grok.shell.warning(`Molecular weight hasn't been calculated for ${errorCount} fragments`);
    });
  }

  detach(): void {
    if ((grok.shell.o as HTMLElement).classList?.contains(MMP_CONTEXT_PANE_CLASS))
      grok.shell.o = ui.div();
    this.lastOpenedHint?.remove();
    super.detach();
  }
}
