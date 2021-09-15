import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { describe } from './describe';

export class SARViewer extends DG.JsViewer {
  initialized: boolean;
  activityColumnName: string;
  activityScalingMethod: string;
  duplicatesHandingMethod: string;
  constructor() {
    super();

    this.activityColumnName = this.string('activityColumnColumnName', 'Activity');
    this.activityScalingMethod = this.string('activityScalingMethod', 'none', {choices: ['none', 'ln', '-ln']});
    this.duplicatesHandingMethod = this.string('duplicatesHandlingMethod', 'median', {choices: ['median']});

    this.initialized = false;
  }

  init() {
    this.initialized = true;
  }

  onTableAttached() {
    this.init();
    
    if (typeof this.dataFrame !== 'undefined') {
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
      this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
      // this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));
    }

    this.render();
  }

  async render() {
    $(this.root).empty();
    if (typeof this.dataFrame !== 'undefined') {
      const grid = await describe(this.dataFrame, this.activityColumnName);
      if (grid !== null) {
        this.root.appendChild(grid.root);
      }
    }
  }
}
