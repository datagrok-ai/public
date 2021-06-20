/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from "cash-dom";

/** A widget that contains the specified HTML code that could be parameterized.
 * All occurrences of ${p1}, ${p2} and ${p3} are replaced with the values of the corresponding properties.
 * */
export class HtmlWidget extends DG.Widget {
  html: string;
  p1: string;
  p2: string;
  p3: string;

  constructor() {
    super(ui.box());

    // properties
    this.html = super.addProperty('html', DG.TYPE.STRING, '<div>P1: ${p1}</div>', {editor: 'markdown'});
    this.p1 = super.addProperty('p1', DG.TYPE.STRING, 'first');
    this.p2 = super.addProperty('p2', DG.TYPE.STRING, '');
    this.p3 = super.addProperty('p3', DG.TYPE.STRING, '');

    this.refresh();
  }

  onPropertyChanged(property: DG.Property | null) {
    this.refresh();
  }

  refresh(): void {
    ui.empty(this.root);

    let s = this.html;
    s = DG.Utils.replaceAll(s, '${p1}', this.p1);
    s = DG.Utils.replaceAll(s, '${p2}', this.p2);
    s = DG.Utils.replaceAll(s, '${p3}', this.p3);

    this.root.appendChild(ui.tools.createElementFromHtml(s));
  }
}

