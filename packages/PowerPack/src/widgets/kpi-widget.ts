/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Property} from "datagrok-api/dg";

export class KpiWidget extends DG.Widget {

  caption: string;
  format: string;
  value: number;

  valueLabel: HTMLDivElement;
  captionLabel: HTMLDivElement;

  constructor() {
    super(ui.div([], 'pp-kpi-host'));

    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Indicator');
    this.value = super.addProperty('value', DG.TYPE.INT, 9999);
    this.format = super.addProperty('format', DG.TYPE.STRING, '###,###,###');

    this.captionLabel = ui.divText('', 'pp-kpi-caption');
    this.valueLabel = ui.divText('', 'pp-kpi-value');
    this.refresh();

    this.root.appendChild(this.captionLabel);
    this.root.appendChild(this.valueLabel);
  }

  refresh(): void {
    let s = DG.format(this.value, this.format);
    this.captionLabel.innerText = this.caption ?? '';
    this.valueLabel.innerText = s;
  }

  onPropertyChanged(property: Property | null): void {
    this.refresh();
  }
}