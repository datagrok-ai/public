import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

function grokViewer(options?: {
  name?: string,
  description?: string,
  icon?: string,
  toolbox?: boolean,
}) {
  return function(constructor: Function) {
  };
}


@grokViewer({
  icon: 'images/1.png',
  toolbox: true,
})
export class FirstViewer extends DG.JsViewer {
  color: string;
  size: number;
  sourceColumnName: string;
  targetColumnName: string;

  constructor() {
    super();
    this.color = this.string('color', 'color');
    this.size = this.int('size', 5);
    this.sourceColumnName = this.string('sourceColumnName', null);
    this.targetColumnName = this.string('targetColumnName', null);
  }

  onPropertyChanged(prop: any) {
    this.render();
  }

  onTableAttached() {
    this.subs.push(this.dataFrame!.selection.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame!.filter.onChanged.subscribe((_) => this.render()));

    this.render();
  }

  render() {
    this.root.innerHTML =
    `First Viewer<br>
    Data properties:<br>
    sourceColumnName: ${this.sourceColumnName}<br>
    targetColumnName: ${this.targetColumnName}<br>
    Style properties:<br>
    color: ${this.color}<br>
    size: ${this.size}<br>
      `;
  }
}
