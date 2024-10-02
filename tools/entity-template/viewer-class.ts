import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

@grok.decorators.viewer({
  name: '#{NAME}',
  description: 'Creates #{NAME} viewer'
})
export class #{NAME} extends DG.JsViewer {
  constructor() {
    super();
  }

  // Additional chart settings
  init() {
  }
  
  // Override to handle property changes
  onPropertyChanged(property : DG.Property | null) {
    super.onPropertyChanged(property);
    this.render();
  }

  // Cancel subscriptions when the viewer is detached
  detach() {
  }
  
  // Stream subscriptions
  onTableAttached() {
    this.subs.push(this.dataFrame!.selection.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame!.filter.onChanged.subscribe((_) => this.render()));

    this.render();
  }

  render() {
    this.root.innerHTML =
      `${this.dataFrame!.toString()}<br>
            Selected: ${this.dataFrame!.selection.trueCount}<br>
            Filtered: ${this.dataFrame!.filter.trueCount}`;
  }
}
