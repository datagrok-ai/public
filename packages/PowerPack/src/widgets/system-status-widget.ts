/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from "datagrok-api/grok";
import {divText} from "datagrok-api/ui";

export class SystemStatusWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.div());

    grok.dapi.admin.getServiceInfos().then((services: any[]) => {
      // let table = ui.table(services, (item, idx) =>
      //   [`${item.key}:`, ui.span([item.status], `grok-plugin-status-${item.status.toLowerCase()}`)])
      // this.root.appendChild(table);

      this.root.appendChild(ui.list(services.map(
        (s) => divText(s.key, `grok-plugin-status-${s.status.toLowerCase()}`))));
    });

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'System');
  }
}

