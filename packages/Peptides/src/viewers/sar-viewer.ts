import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {PeptidesController} from '../peptides';
import {SARMultipleFilter} from '../utils/SAR-multiple-filter';
// import {PeptidesModel} from '../model';

/**
 * Structure-activity relationship viewer.
 */
export class SARViewer extends DG.JsViewer {
  protected viewerGrid: DG.Grid | null;
  protected sourceGrid: DG.Grid | null;
  protected activityColumnName: string;
  protected scaling: string;
  protected bidirectionalAnalysis: boolean;
  protected filterMode: boolean;
  protected statsDf: DG.DataFrame | null;
  protected initialized: boolean;
  protected viewGridInitialized: boolean;
  protected aminoAcidResidue;
  protected viewerVGrid: DG.Grid | null;
  protected grouping: boolean;
  protected groupMapping: StringDictionary | null;
  // model: PeptidesModel | null;
  protected _name: string = 'Monomer-Positions';
  protected controller: PeptidesController | null;
  protected multipleFilter: SARMultipleFilter;
  // protected pValueThreshold: number;
  // protected amountOfBestAARs: number;
  // duplicatesHandingMethod: string;

  constructor() {
    super();

    this.viewerGrid = null;
    this.viewerVGrid = null;
    this.statsDf = null;
    this.groupMapping = null;
    this.initialized = false;
    this.aminoAcidResidue = 'AAR';
    this.viewGridInitialized = false;
    // this.model = null;
    this.controller = null;

    //TODO: find a way to restrict activityColumnName to accept only numerical columns (double even better)
    this.activityColumnName = this.string('activityColumnName');
    this.scaling = this.string('scaling', 'none', {choices: ['none', 'lg', '-lg']});
    this.filterMode = this.bool('filterMode', false);
    this.bidirectionalAnalysis = this.bool('bidirectionalAnalysis', false);
    this.grouping = this.bool('grouping', false);

    this.sourceGrid = null;
    this.multipleFilter = new SARMultipleFilter(this.filterMode);
    this.multipleFilter.addResource('residueColumnName', this.aminoAcidResidue);
  }

  get name() {
    return this._name;
  }

  /*private get mask() {
    const df = this.dataFrame!;
    return this.filterMode ? df.filter : df.selection;
  }

  private set mask(mask: DG.BitSet) {
    const df = this.dataFrame!;
    if (this.filterMode)
      df.filter = mask;
    else
      df.rows.select((row) => mask.get(row.idx));
  }*/

  private get filter() {
    return this.dataFrame!.filter;// ?? DG.BitSet.create(1, (_) => true);
  }

  init() {
  }

  async onTableAttached() {
    this.sourceGrid = this.view?.grid ?? (grok.shell.v as DG.TableView).grid;
    this.dataFrame?.setTag('dataType', 'peptides');
    this.controller = await PeptidesController.getInstance(this.dataFrame!);
    this.initialized = true;
    // this.model = PeptidesModel.getOrInit(this.dataFrame!);
    // this.model = this.controller.getOrInitModel();

    this.multipleFilter.addResource('dataFrame', this.dataFrame);
    this.multipleFilter.resetSelection();

    this.subs.push(this.controller.onStatsDataFrameChanged.subscribe((data) => this.statsDf = data));
    this.subs.push(this.controller.onSARGridChanged.subscribe((data) => {
      this.viewerGrid = data;
      this.multipleFilter.addResource('grid', this.viewerGrid);
      this.multipleFilter.onSARGridChanged();
      this.render(false);
    }));
    this.subs.push(this.controller.onSARVGridChanged.subscribe((data) => this.viewerVGrid = data));
    this.subs.push(this.controller.onGroupMappingChanged.subscribe((data) => {
      this.groupMapping = data;
      this.multipleFilter.addResource('groupMapping', this.groupMapping);
    }));
    this.subs.push(this.dataFrame?.onRowsFiltering.subscribe((_) => {
      this.multipleFilter.maskRows();
    })!);
    this.subs.push(this.dataFrame!.onSelectionChanged.subscribe((_) => {
      this.multipleFilter.maskRows();
    }));

    await this.render();
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /**
   * Function that is executed when the property is changed.
   *
   * @param {DG.Property} property New property.
   * @memberof SARViewer
   */
  async onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);

    if (!this.initialized) {
      this.init();
      return;
    }

    if (property.name === 'grouping')
      this.multipleFilter.resetSelection();


    if (property.name === 'filterMode')
      this.multipleFilter.filteringMode = this.filterMode;

    if (property.name === 'activityColumnName')
      this.multipleFilter.addResource('activityColumnName', `${this.activityColumnName}Scaled`);

    if (property.name === 'scaling' && typeof this.dataFrame !== 'undefined') {
      const minActivity = DG.Stats.fromColumn(
        this.dataFrame!.col(this.activityColumnName)!,
        this.dataFrame.filter,
      ).min;
      if (minActivity && minActivity <= 0 && this.scaling !== 'none') {
        grok.shell.warning(`Could not apply ${this.scaling}: ` +
          `activity column ${this.activityColumnName} contains zero or negative values, falling back to 'none'.`);
        property.set(this, 'none');
        return;
      }
    }

    await this.render();
  }

  applyBitset() {
    if (!this.viewerGrid)
      return;

    const viewerGridDf = this.viewerGrid!.dataFrame;

    if (viewerGridDf && viewerGridDf.currentCell.value && viewerGridDf.currentCol.name !== this.aminoAcidResidue) {
      //const currentAAR: string = viewerGridDf.get(this.aminoAcidResidue, viewerGridDf.currentRowIdx);
      //const currentPosition = viewerGridDf.currentCol.name;
      //const aarLabel = `${currentAAR === '-' ? 'Gap' : currentAAR} - ${currentPosition}`;
      const splitColName = '~splitCol';
      const otherLabel = 'Other';
      const aarLabel = this.multipleFilter.filterLabel;
      const splitCol = this.dataFrame!.col(splitColName) ?? this.dataFrame!.columns.addNew(splitColName, 'string');

      (splitCol! as DG.Column).init((i) => this.filter.get(i) ? aarLabel : otherLabel);

      const colorMap: {[index: string]: string | number} = {};

      colorMap[otherLabel] = DG.Color.blue;
      colorMap[aarLabel] = DG.Color.orange;
      // colorMap[currentAAR] = cp.getColor(currentAAR);
      this.dataFrame!.getCol(splitColName).colors.setCategorical(colorMap);
    }
  }

  accordionFunc(accordion: DG.Accordion) {
    if (accordion.context instanceof DG.RowGroup) {
      const originalDf: DG.DataFrame = DG.toJs(accordion.context.dataFrame);
      const viewerDf = this.viewerGrid!.dataFrame;

      if (
        originalDf.getTag('dataType') === 'peptides' &&
        originalDf.col('~splitCol') &&
        viewerDf &&
        viewerDf.currentCol !== null
      ) {
        const labelStr = this.multipleFilter.filterLabel;
        const currentColor = DG.Color.toHtml(DG.Color.orange);
        const otherColor = DG.Color.toHtml(DG.Color.blue);
        const currentLabel = ui.label(labelStr, {style: {color: currentColor}});
        const otherLabel = ui.label('Other', {style: {color: otherColor}});
        const elements: (HTMLLabelElement | HTMLElement)[] = [currentLabel, otherLabel];
        const distPane = accordion.getPane('Distribution');

        if (distPane)
          accordion.removePane(distPane);

        const getContent = () => {
          const hist = originalDf.clone(this.dataFrame!.filter).plot.histogram({
          // const hist = originalDf.plot.histogram({
            filteringEnabled: false,
            valueColumnName: `${this.activityColumnName}Scaled`,
            splitColumnName: '~splitCol',
            legendVisibility: 'Never',
            showXAxis: true,
            showColumnSelector: false,
            showRangeSlider: false,
          }).root;

          hist.style.width = 'auto';
          elements.push(hist);

          this.multipleFilter.addResource('activityColumnName', `${this.activityColumnName}Scaled`);

          const stats = this.multipleFilter.getStatistics();
          const tableMap: StringDictionary = {
            'Statistics:': '',
            'Count': stats.count.toString(),
            'p-value': stats.pValue < 0.01 ? '<0.01' : stats.pValue.toFixed(2),
            'Mean difference': stats.meanDifference.toFixed(2),
          };

          elements.push(ui.tableFromMap(tableMap));
          return ui.divV(elements);
        };

        accordion.addPane('Distribution', getContent.bind(this), true);
      }
    }
  }

  /**
   * Viewer render function.
   *
   * @param {boolean} [computeData=true] Recalculate data.
   * @memberof SARViewer
   */
  async render(computeData = true) {
    if (!this.initialized)
      return;

    try {
    //TODO: optimize. Don't calculate everything again if only view changes
      if (typeof this.dataFrame !== 'undefined' && this.activityColumnName && this.sourceGrid) {
        if (computeData) {
          await this.controller!.updateData(this.activityColumnName, this.scaling, this.sourceGrid,
            this.bidirectionalAnalysis, this.filter, this.grouping);
        }

        if (this.viewerGrid !== null && this.viewerVGrid !== null) {
          $(this.root).empty();
          const title = ui.h1(this._name, {style: {'align-self': 'center'}});
          const gridRoot = this.viewerGrid.root;
          gridRoot.style.width = 'auto';
          this.root.appendChild(ui.divV([title, gridRoot]));
        this.viewerGrid.dataFrame!.onCurrentCellChanged.subscribe((_) => {
          this.applyBitset();
          syncGridsFunc(false, this.viewerGrid!, this.viewerVGrid!, this.aminoAcidResidue);
        });
        this.viewerVGrid.dataFrame!.onCurrentCellChanged.subscribe((_) => {
          syncGridsFunc(true, this.viewerGrid!, this.viewerVGrid!, this.aminoAcidResidue);
        });
        grok.events.onAccordionConstructed.subscribe((accordion: DG.Accordion) => {
          this.accordionFunc(accordion);
        });
        }
      }
      //fixes viewers not rendering immediately after analyze.
      this.viewerGrid?.invalidate();
    } catch (error) {
      console.warn(error);
    }
  }
}

/**
 * Vertical structure activity relationship viewer.
 */
export class SARViewerVertical extends DG.JsViewer {
  viewerVGrid: DG.Grid | null;
  // model: PeptidesModel | null;
  protected _name = 'Sequence-Activity relationship';
  controller: PeptidesController | null;

  constructor() {
    super();

    this.viewerVGrid = null;
    this.controller = null;
  }

  get name() {
    return this._name;
  }

  async onTableAttached() {
    // this.model = PeptidesModel.getOrInit(this.dataFrame!);
    this.controller = await PeptidesController.getInstance(this.dataFrame!);

    this.subs.push(this.controller.onSARVGridChanged.subscribe((data) => {
      this.viewerVGrid = data;
      this.render();
    }));
  }

  render() {
    if (this.viewerVGrid) {
      $(this.root).empty();
      this.root.appendChild(this.viewerVGrid.root);
    }
    this.viewerVGrid?.invalidate();
  }
}

//TODO: refactor, move
function syncGridsFunc(sourceVertical: boolean, viewerGrid: DG.Grid, viewerVGrid: DG.Grid, aminoAcidResidue: string) {
  if (viewerGrid && viewerGrid.dataFrame && viewerVGrid && viewerVGrid.dataFrame) {
    if (sourceVertical) {
      const dfCell = viewerVGrid.dataFrame.currentCell;
      if (dfCell.column === null || dfCell.column.name !== 'Diff')
        return;

      const otherColName: string = viewerVGrid.dataFrame.get('Pos', dfCell.rowIndex);
      const otherRowName: string = viewerVGrid.dataFrame.get(aminoAcidResidue, dfCell.rowIndex);
      let otherRowIndex = -1;
      for (let i = 0; i < viewerGrid.dataFrame.rowCount; i++) {
        if (viewerGrid.dataFrame.get(aminoAcidResidue, i) === otherRowName) {
          otherRowIndex = i;
          break;
        }
      }
      if (otherRowIndex !== -1)
        viewerGrid.dataFrame.currentCell = viewerGrid.dataFrame.cell(otherRowIndex, otherColName);
    } else {
      const otherPos: string = viewerGrid.dataFrame.currentCol?.name;
      if (typeof otherPos === 'undefined' && otherPos !== aminoAcidResidue)
        return;

      const otherAAR: string =
        viewerGrid.dataFrame.get(aminoAcidResidue, viewerGrid.dataFrame.currentRowIdx);
      let otherRowIndex = -1;
      for (let i = 0; i < viewerVGrid.dataFrame.rowCount; i++) {
        if (
          viewerVGrid.dataFrame.get(aminoAcidResidue, i) === otherAAR &&
          viewerVGrid.dataFrame.get('Pos', i) === otherPos
        ) {
          otherRowIndex = i;
          break;
        }
      }
      if (otherRowIndex !== -1)
        viewerVGrid.dataFrame.currentCell = viewerVGrid.dataFrame.cell(otherRowIndex, 'Diff');
    }
  }
}
