import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {describe} from './describe';
import $ from "cash-dom";

export class SARViewer extends DG.JsViewer {
  private grid: DG.Viewer | null;
  protected activityColumnColumnName: string;
  protected activityScalingMethod: string;
  protected filterMode: boolean;
  // protected showHistogram: boolean;
  // duplicatesHandingMethod: string;
  constructor() {
    super();

    this.grid = null;

    //TODO: find a way to restrict activityColumnColumnName to accept only numerical columns (double even better)
    this.activityColumnColumnName = this.string('activityColumnColumnName');
    this.activityScalingMethod = this.string('activityScalingMethod', 'none', {choices: ['none', 'lg', '-lg']});
    this.filterMode = this.bool('filterMode', false);
    // this.showHistogram = this.bool('showHistogram', false);
    // this.duplicatesHandingMethod = this.string('duplicatesHandlingMethod', 'median', {choices: ['median']});
  }

  async onTableAttached() {
    if (typeof this.dataFrame !== 'undefined') {
      // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_: any) => this.render()));
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
      this.grid = await describe(this.dataFrame, this.activityColumnColumnName, this.activityScalingMethod, this.filterMode);

      if (this.grid !== null) {
        $(this.root).empty();
        this.root.appendChild(this.grid.root);

        // if (this.showHistogram) {
        //   this.root.appendChild(grid.table.plot.histogram().root);
        // }
      }
    }
  }
}
