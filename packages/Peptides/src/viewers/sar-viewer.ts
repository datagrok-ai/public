import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {PeptidesController} from '../peptides';
import {SARMultipleFilter} from '../utils/SAR-multiple-filter';
import * as C from '../utils/constants';
// import {PeptidesModel} from '../model';

/**
 * Structure-activity relationship viewer.
 */
export class SARViewer extends DG.JsViewer {
  viewerGrid: DG.Grid | null;
  sourceGrid: DG.Grid | null;
  // activityColumnName: string;
  scaling: string;
  bidirectionalAnalysis: boolean;
  filterMode: boolean;
  statsDf: DG.DataFrame | null;
  initialized: boolean;
  viewGridInitialized: boolean;
  viewerVGrid: DG.Grid | null;
  grouping: boolean;
  groupMapping: StringDictionary | null;
  _titleHost = ui.divText('Monomer-Positions', {id: 'pep-viewer-title'});
  _name = 'Structure-Activity Relationship';
  controller: PeptidesController | null;
  multipleFilter: SARMultipleFilter;
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
    this.viewGridInitialized = false;
    // this.model = null;
    this.controller = null;

    //TODO: find a way to restrict activityColumnName to accept only numerical columns (double even better)
    // this.activityColumnName = this.string('activityColumnName');
    this.scaling = this.string('scaling', 'none', {choices: ['none', 'lg', '-lg']});
    this.filterMode = this.bool('filterMode', false);
    this.bidirectionalAnalysis = this.bool('bidirectionalAnalysis', false);
    this.grouping = this.bool('grouping', false);

    this.sourceGrid = null;
    this.multipleFilter = new SARMultipleFilter(this.filterMode);
    this.multipleFilter.addResource('residueColumnName', C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
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
  get selection() {
    return this.dataFrame!.selection;
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
    this.multipleFilter.addResource('activityColumnName', C.COLUMNS_NAMES.ACTIVITY_SCALED);

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

    if (!this.initialized)
      return;

    if (property.name === 'grouping')
      this.multipleFilter.resetSelection();


    if (property.name === 'filterMode')
      this.multipleFilter.filteringMode = this.filterMode;

    if (property.name === 'activityColumnName')
      this.multipleFilter.addResource('activityColumnName', C.COLUMNS_NAMES.ACTIVITY_SCALED);

    if (property.name === 'scaling' && typeof this.dataFrame !== 'undefined') {
      const minActivity = DG.Stats.fromColumn(this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY), this.dataFrame.filter)
        .min;
      if (minActivity && minActivity <= 0 && this.scaling !== 'none') {
        grok.shell.warning(`Could not apply ${this.scaling}: ` +
          `activity column ${C.COLUMNS_NAMES.ACTIVITY} contains zero or negative values, falling back to 'none'.`);
        property.set(this, 'none');
        return;
      }
    }

    await this.render();
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


    //TODO: optimize. Don't calculate everything again if only view changes
    if (typeof this.dataFrame !== 'undefined' && this.sourceGrid) {
      if (computeData) {
        await this.controller!.updateData(this.scaling, this.sourceGrid, this.bidirectionalAnalysis, this.filter,
          this.grouping);
      }

      if (this.viewerGrid !== null && this.viewerVGrid !== null) {
        $(this.root).empty();
        const gridRoot = this.viewerGrid.root;
        gridRoot.style.width = 'auto';
        this.root.appendChild(ui.divV([this._titleHost, gridRoot]));
      this.viewerGrid.dataFrame!.onCurrentCellChanged.subscribe((_) => {
        // this.applyBitset();
        syncGridsFunc(false, this.viewerGrid!, this.viewerVGrid!);
      });
      this.viewerVGrid.dataFrame!.onCurrentCellChanged.subscribe((_) => {
        syncGridsFunc(true, this.viewerGrid!, this.viewerVGrid!);
      });
      // grok.events.onAccordionConstructed.subscribe((accordion: DG.Accordion) => {
      //   this.accordionFunc(accordion);
      // });
      }
    }
    //fixes viewers not rendering immediately after analyze.
    this.viewerGrid?.invalidate();
  }
}

/**
 * Vertical structure activity relationship viewer.
 */
export class SARViewerVertical extends DG.JsViewer {
  viewerVGrid: DG.Grid | null;
  // model: PeptidesModel | null;
  protected _name = 'Sequence-Activity relationship';
  _titleHost = ui.divText('Most Potent Residues', {id: 'pep-viewer-title'});
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
      this.root.appendChild(ui.divV([this._titleHost, this.viewerVGrid.root]));
    }
    this.viewerVGrid?.invalidate();
  }
}

//TODO: refactor, move
function syncGridsFunc(sourceVertical: boolean, viewerGrid: DG.Grid, viewerVGrid: DG.Grid) {
  if (viewerGrid && viewerGrid.dataFrame && viewerVGrid && viewerVGrid.dataFrame) {
    if (sourceVertical) {
      const dfCell = viewerVGrid.dataFrame.currentCell;
      if (dfCell.column === null || dfCell.column.name !== 'Diff')
        return;

      const otherColName: string = viewerVGrid.dataFrame.get('Pos', dfCell.rowIndex);
      const otherRowName: string = viewerVGrid.dataFrame.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, dfCell.rowIndex);
      let otherRowIndex = -1;
      for (let i = 0; i < viewerGrid.dataFrame.rowCount; i++) {
        if (viewerGrid.dataFrame.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, i) === otherRowName) {
          otherRowIndex = i;
          break;
        }
      }
      if (otherRowIndex !== -1)
        viewerGrid.dataFrame.currentCell = viewerGrid.dataFrame.cell(otherRowIndex, otherColName);
    } else {
      const otherPos: string = viewerGrid.dataFrame.currentCol?.name;
      if (typeof otherPos === 'undefined' && otherPos !== C.COLUMNS_NAMES.AMINO_ACID_RESIDUE)
        return;

      const otherAAR: string =
        viewerGrid.dataFrame.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, viewerGrid.dataFrame.currentRowIdx);
      let otherRowIndex = -1;
      for (let i = 0; i < viewerVGrid.dataFrame.rowCount; i++) {
        if (
          viewerVGrid.dataFrame.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, i) === otherAAR &&
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
