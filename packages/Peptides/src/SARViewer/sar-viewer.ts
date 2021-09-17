import * as DG from 'datagrok-api/dg';
import {describe} from './describe';

export class SARViewer extends DG.JsViewer {
  protected activityColumnColumnName: string;
  protected activityScalingMethod: string;
  protected showHistogram: boolean;
  // duplicatesHandingMethod: string;
  constructor() {
    super();

    //TODO: find a way to restrict activityColumnColumnName to accept only numerical columns (double even better)
    this.activityColumnColumnName = this.string('activityColumnColumnName');
    this.activityScalingMethod = this.string('activityScalingMethod', 'none', {choices: ['none', 'lg', '-lg']});
    this.showHistogram = this.bool('showHistogram', false);
    // this.duplicatesHandingMethod = this.string('duplicatesHandlingMethod', 'median', {choices: ['median']});
  }

  async onTableAttached() {
    if (typeof this.dataFrame !== 'undefined') {
      // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_: any) => this.render()));
      // this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_: any) => this.render()));
    }

    this.render();
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);
    this.render();
  }

  async render() {
    if (typeof this.dataFrame !== 'undefined' && this.activityColumnColumnName !== null) {
      const grid = await describe(this.dataFrame, this.activityColumnColumnName, this.activityScalingMethod);

      if (grid !== null) {
        $(this.root).empty();
        this.root.appendChild(grid.root);

        if (this.showHistogram) {
          this.root.appendChild(grid.table.plot.histogram().root);
        }
      }
    }
  }
}
