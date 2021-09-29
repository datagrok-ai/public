import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {describe} from './describe';
import $ from 'cash-dom';
import {tTest} from '../utils/misc';

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
    const accordionFunc = (accordion: DG.Accordion) => {
      if (accordion.context instanceof DG.DataFrame || typeof accordion.context.dataFrame !== 'undefined') {
        // console.log(accordion.context);
        const originalDf: DG.DataFrame = accordion.context instanceof DG.DataFrame ? accordion.context : accordion.context.dataFrame;

        if (originalDf.getTag('dataType') === 'peptides' && originalDf.col('~splitCol')) {
          let histPane = accordion.getPane('Distribution');
          histPane = histPane ? histPane : accordion.addPane('Distribution', () => {
            return originalDf.plot.histogram({
              value: `~${this.activityColumnColumnName}Scaled`,
              'splitColumnName': '~splitCol',
            }).root;
          }, true);
          // console.log('from accordion');
          // console.log(originalDf.columns.names());
          // console.log(originalDf.col('~splitCol')?.categories);
        }
      }
    };
    // grok.events.onAccordionConstructed.subscribe();
    if (typeof this.dataFrame !== 'undefined') {
    //   this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_: any) => this.render()));
    //   this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_: any) => this.render()));
    //   this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(grok.events.onAccordionConstructed, 50).subscribe(accordionFunc));
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

      // const a1 = [13.3, 6.0, 20.0, 8.0, 14.0, 19.0, 18.0, 25.0, 16.0, 24.0, 15.0, 1.0, 15.0];
      // const a2 = [22.0, 16.0, 21.7, 21.0, 30.0, 26.0, 12.0, 23.2, 28.0, 23.0];
      // console.log(tTest(a1, a2, 0.05, false, false));
      // console.log(tTest(a1, a2, 0.05, true, false));
      // console.log(tTest(a1, a2, 0.05, true, true));
    }
  }
}
