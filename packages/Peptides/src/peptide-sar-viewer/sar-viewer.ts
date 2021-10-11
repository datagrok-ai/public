import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {describe} from './describe';

export class SARViewer extends DG.JsViewer {
  private grid: DG.Viewer | null;
  protected activityColumnColumnName: string;
  protected activityScalingMethod: string;
  // protected filterMode: boolean;
  protected statsDf: DG.DataFrame | null;
  protected initialized: boolean;
  // duplicatesHandingMethod: string;
  constructor() {
    super();

    this.grid = null;
    this.statsDf = null;
    this.initialized = false;

    //TODO: find a way to restrict activityColumnColumnName to accept only numerical columns (double even better)
    this.activityColumnColumnName = this.string('activityColumnColumnName');
    this.activityScalingMethod = this.string('activityScalingMethod', 'none', {choices: ['none', 'lg', '-lg']});
    // this.filterMode = this.bool('filterMode', false);
    // this.duplicatesHandingMethod = this.string('duplicatesHandlingMethod', 'median', {choices: ['median']});
  }

  init() {
    this.initialized = true;
  }

  onTableAttached() {
    const accordionFunc = (accordion: DG.Accordion) => {
      if (accordion.context instanceof DG.RowGroup) {
        const originalDf: DG.DataFrame = accordion.context.dataFrame;

        if (originalDf.getTag('dataType') === 'peptides' && originalDf.col('~splitCol')) {
          const currentAAR: string = this.grid?.table.get('aminoAcidResidue', this.grid?.table.currentRowIdx);
          const currentPosition = this.grid?.table.currentCol.name;

          const labelStr = `${currentAAR === '-' ? 'Empty' : currentAAR} - ${currentPosition}`;
          // const currentColor = DG.Color.toHtml(DG.Color.getCategoryColor(originalDf.getCol('~splitCol'), labelStr));
          // const otherColor = DG.Color.toHtml(DG.Color.getCategoryColor(originalDf.getCol('~splitCol'), 'Other'));
          const currentLabel = ui.label(labelStr, {style: {color: DG.Color.toHtml(DG.Color.orange)}});
          const otherLabel = ui.label('Other', {style: {color: DG.Color.toHtml(DG.Color.blue)}});

          const elements: (HTMLLabelElement | HTMLElement)[] = [currentLabel, otherLabel];

          accordion.getPane('Distribution') ? null : accordion.addPane('Distribution', () => {
            const hist = originalDf.plot.histogram({
              valueColumnName: `~${this.activityColumnColumnName}Scaled`,
              splitColumnName: '~splitCol',
              legendVisibility: 'Never',
              showXAxis: true,
              showColumnSelector: false,
              showRangeSlider: false,
            }).root;
            elements.push(hist);

            const tableMap: {[key: string]: string} = {'Statistics:': ''};
            for (const colName of new Set(['Count', 'p-value', 'Mean difference'])) {
              const query = `aminoAcidResidue = ${currentAAR} and position = ${currentPosition}`;
              const text = `${this.statsDf?.groupBy([colName]).where(query).aggregate().get(colName, 0).toFixed(2)}`;
              tableMap[colName] = text;
            }
            elements.push(ui.tableFromMap(tableMap));

            return ui.divV(elements);
          }, true);
        }
      }
    };

    this.subs.push(DG.debounce(grok.events.onAccordionConstructed, 50).subscribe(accordionFunc));

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

    if (property.name === 'activityScalingMethod') {
      const minActivity = this.dataFrame?.col(this.activityColumnColumnName)?.min;
      if (minActivity && minActivity <= 0 && this.activityScalingMethod !== 'none') {
        grok.shell.warning(`Could not apply ${this.activityScalingMethod}: ` +
          `activity column ${this.activityColumnColumnName} contains zero or negative values, falling back to 'none'.`);
        property.set(this, 'none');
        return;
      }
    }

    this.render();
  }

  async render() {
    if (!this.initialized) {
      return;
    }
    //TODO: optimize. Don't calculate everything again if only view changes
    if (typeof this.dataFrame !== 'undefined' && this.activityColumnColumnName) {
      [this.grid, this.statsDf] = await describe(
        this.dataFrame,
        this.activityColumnColumnName,
        this.activityScalingMethod,
        // this.filterMode,
        false,
      );

      if (this.grid !== null) {
        $(this.root).empty();
        this.root.appendChild(this.grid.root);
      }
    }
  }
}
