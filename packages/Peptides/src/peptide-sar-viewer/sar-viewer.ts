import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {describe} from './describe';
import $ from 'cash-dom';

export class SARViewer extends DG.JsViewer {
  private grid: DG.Viewer | null;
  protected activityColumnColumnName: string;
  protected activityScalingMethod: string;
  protected filterMode: boolean;
  // duplicatesHandingMethod: string;
  constructor() {
    super();

    this.grid = null;

    //TODO: find a way to restrict activityColumnColumnName to accept only numerical columns (double even better)
    this.activityColumnColumnName = this.string('activityColumnColumnName');
    this.activityScalingMethod = this.string('activityScalingMethod', 'none', {choices: ['none', 'lg', '-lg']});
    this.filterMode = this.bool('filterMode', false);
    // this.duplicatesHandingMethod = this.string('duplicatesHandlingMethod', 'median', {choices: ['median']});
  }

  async onTableAttached() {
    if (typeof this.dataFrame !== 'undefined') {
      // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_: any) => this.render()));
      // this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_: any) => this.render()));
    }

    this.render();
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);
    this.render();
  }

  async render() {
    //TODO: optimize. Don't calculate everything again if only view changes
    if (typeof this.dataFrame !== 'undefined' && this.activityColumnColumnName !== null) {
      this.grid = await describe(
        this.dataFrame,
        this.activityColumnColumnName,
        this.activityScalingMethod,
        this.filterMode,
      );

      if (this.grid !== null) {
        $(this.root).empty();
        this.root.appendChild(this.grid.root);
      }

      grok.events.onAccordionConstructed.subscribe((accordion: DG.Accordion) => {
        if (accordion.context instanceof DG.DataFrame || typeof accordion.context.dataFrame !== 'undefined') {
          const originalDf: DG.DataFrame = accordion.context instanceof DG.DataFrame ? accordion.context : DG.toJs(accordion.context.dataFrame);

          if (originalDf.getTag('dataType') === 'peptides') {
            // const currentDf = originalDf.clone(originalBitset);

            // //FIXME: Can this be done better?
            // //add bool col to split df on selection/filter instead
            // let otherDf = originalDf.clone();
            // const otherBitset = this.filterMode ? otherDf.filter : otherDf.selection;
            // otherBitset.init((i) => !originalBitset.get(i));
            // otherDf = otherDf.clone(otherBitset);
    
            let histPane = accordion.getPane('Histogram');
            histPane = histPane ? histPane : accordion.addPane('Histogram', () => {
              return originalDf.plot.histogram({value: this.activityColumnColumnName + 'Scaled', 'splitColumnName': 'splitCol'});
            }, true);
          }
        }
      });

      const columnNames = this.dataFrame.columns.names();
      let tempCol = this.dataFrame.col('splitCol');
      tempCol ? columnNames.splice(columnNames.indexOf(tempCol.name), 1) : null;
      tempCol = this.dataFrame.col(this.activityColumnColumnName + 'Scaled');
      tempCol ? columnNames.splice(columnNames.indexOf(tempCol.name), 1) : null;

      //@ts-ignore: *sigh*
      grok.shell.v.grid.columns.setVisible(columnNames);
    }
  }
}
