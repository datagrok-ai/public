import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {describe} from './describe';

export class SARViewer extends DG.JsViewer {
  private viewerGrid: DG.Grid | null;
  private sourceGrid: DG.Grid;
  private progress: DG.TaskBarProgressIndicator;
  protected activityColumnColumnName: string;
  protected activityScalingMethod: string;
  protected bidirectionalAnalysis: boolean;
  protected filterMode: boolean;
  protected statsDf: DG.DataFrame | null;
  protected initialized: boolean;
  protected viewGridInitialized: boolean;
  protected aminoAcidResidue;
  protected _initialBitset: DG.BitSet | null;
  // duplicatesHandingMethod: string;
  constructor() {
    super();
    this.progress = DG.TaskBarProgressIndicator.create('Loading SAR viewer');

    this.viewerGrid = null;
    this.statsDf = null;
    this.initialized = false;
    this.aminoAcidResidue = 'aminoAcidResidue';
    this._initialBitset = null;
    this.viewGridInitialized = false;

    //TODO: find a way to restrict activityColumnColumnName to accept only numerical columns (double even better)
    this.activityColumnColumnName = this.string('activityColumnColumnName');
    this.activityScalingMethod = this.string('activityScalingMethod', 'none', {choices: ['none', 'lg', '-lg']});
    this.filterMode = this.bool('filterMode', false);
    this.bidirectionalAnalysis = this.bool('bidirectionalAnalysis', false);
    // this.duplicatesHandingMethod = this.string('duplicatesHandlingMethod', 'median', {choices: ['median']});

    this.sourceGrid = (grok.shell.v as DG.TableView).grid;
  }

  init() {
    this._initialBitset = this.dataFrame!.filter.clone();
    this.initialized = true;
  }

  onTableAttached() {
    this.render();
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);

    if (!this.initialized) {
      this.init();
      return;
    }

    if (property.name === 'activityScalingMethod' && typeof this.dataFrame !== 'undefined') {
      const minActivity = DG.Stats.fromColumn(
        this.dataFrame.col(this.activityColumnColumnName)!,
        this._initialBitset,
      ).min;
      if (minActivity && minActivity <= 0 && this.activityScalingMethod !== 'none') {
        grok.shell.warning(`Could not apply ${this.activityScalingMethod}: ` +
          `activity column ${this.activityColumnColumnName} contains zero or negative values, falling back to 'none'.`);
        property.set(this, 'none');
        return;
      }
    }

    this.render();
  }

  private applyBitset() {
    if (
      this.dataFrame &&
      this.viewerGrid &&
      this.viewerGrid.table.currentCell.value &&
      this.viewerGrid.table.currentCol.name !== this.aminoAcidResidue
    ) {
      const currentAAR: string = this.viewerGrid.table.get(this.aminoAcidResidue, this.viewerGrid.table.currentRowIdx);
      const currentPosition = this.viewerGrid.table.currentCol.name;

      const splitColName = '~splitCol';
      const otherLabel = 'Other';
      const aarLabel = `${currentAAR === '-' ? 'Empty' : currentAAR} - ${currentPosition}`;

      if (!this.dataFrame.col(splitColName)) {
        this.dataFrame.columns.addNew(splitColName, 'string');
      }

      const isChosen = (i: number) => this.dataFrame!.get(currentPosition, i) === currentAAR;
      this.dataFrame.getCol(splitColName).init((i) => isChosen(i) ? aarLabel : otherLabel);

      if (this.filterMode) {
        this.dataFrame.selection.setAll(false, false);
        this.dataFrame.filter.init(isChosen).and(this._initialBitset!, false);
      } else {
        this.dataFrame.filter.copyFrom(this._initialBitset!);
        this.dataFrame.selection.init(isChosen).and(this._initialBitset!, false);
      }

      // df.getCol(splitColName).setCategoryOrder([otherLabel, aarLabel]);
      const colorMap: {[index: string]: string | number} = {};
      colorMap[otherLabel] = DG.Color.blue;
      colorMap[aarLabel] = DG.Color.orange;
      // colorMap[currentAAR] = cp.getColor(currentAAR);
      this.dataFrame.getCol(splitColName).colors.setCategorical(colorMap);
    }
  }

  private accordionFunc(accordion: DG.Accordion) {
    if (accordion.context instanceof DG.RowGroup) {
      const originalDf: DG.DataFrame = accordion.context.dataFrame;

      if (originalDf.getTag('dataType') === 'peptides' && originalDf.col('~splitCol')) {
        const currentAAR: string = this.viewerGrid?.table.get(
          this.aminoAcidResidue,
          this.viewerGrid?.table.currentRowIdx,
        );
        const currentPosition = this.viewerGrid?.table.currentCol.name;

        const labelStr = `${currentAAR === '-' ? 'Empty' : currentAAR} - ${currentPosition}`;
        const currentColor = DG.Color.toHtml(DG.Color.orange);
        const otherColor = DG.Color.toHtml(DG.Color.blue);
        const currentLabel = ui.label(labelStr, {style: {color: currentColor}});
        const otherLabel = ui.label('Other', {style: {color: otherColor}});

        const elements: (HTMLLabelElement | HTMLElement)[] = [currentLabel, otherLabel];

        const distPane = accordion.getPane('Distribution');
        if (distPane) {
          accordion.removePane(distPane);
        }
        accordion.addPane('Distribution', () => {
          const hist = originalDf.clone(this._initialBitset).plot.histogram({
            valueColumnName: `${this.activityColumnColumnName}Scaled`,
            splitColumnName: '~splitCol',
            legendVisibility: 'Never',
            showXAxis: true,
            showColumnSelector: false,
            showRangeSlider: false,
          }).root;
          elements.push(hist);

          const tableMap: {[key: string]: string} = {'Statistics:': ''};
          for (const colName of new Set(['Count', 'p-value', 'Mean difference'])) {
            const query = `${this.aminoAcidResidue} = ${currentAAR} and position = ${currentPosition}`;
            const textNum = this.statsDf?.groupBy([colName]).where(query).aggregate().get(colName, 0);
            const text = textNum === 0 ? '<0.01' : `${colName === 'Count' ? textNum : textNum.toFixed(2)}`;
            tableMap[colName] = text;
          }
          elements.push(ui.tableFromMap(tableMap));

          return ui.divV(elements);
        }, true);
      }
    }
  };

  async render() {
    if (!this.initialized) {
      return;
    }
    //TODO: optimize. Don't calculate everything again if only view changes
    if (typeof this.dataFrame !== 'undefined' && this.activityColumnColumnName) {
      [this.viewerGrid, this.statsDf] = await describe(
        this.dataFrame,
        this.activityColumnColumnName,
        this.activityScalingMethod,
        this.sourceGrid,
        this.bidirectionalAnalysis,
        this._initialBitset,
      );

      if (this.viewerGrid !== null) {
        $(this.root).empty();
        this.root.appendChild(this.viewerGrid.root);
        this.viewerGrid.table.onCurrentCellChanged.subscribe((_) => this.applyBitset());

        grok.events.onAccordionConstructed.subscribe((accordion: DG.Accordion) => this.accordionFunc(accordion));
      }
    }
    this.progress.close();
  }
}
