import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class DGMarkdown extends HTMLElement {
  private _content: string = '';

  constructor() {
    super();
  }

  public set content(val: string) {
    this._content = val;
    this.render();
  }

  private render() {
    ui.empty(this);
    const markdown = ui.markdown(this._content);
    this.appendChild(markdown);
  }
}
