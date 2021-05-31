/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class ServicesWidget extends DG.Widget {
  caption: string;

  constructor() {
    let root = ui.div();
    super(root);

    grok.dapi.admin.getServiceInfos().then((services: any[]) => {
      console.log(services);
      let table = ui.table(services, (item, idx) =>
        [`${item.key}:`, ui.span([item.status], `grok-plugin-status-${item.status.toLowerCase()}`)])
      root.appendChild(table);
    });

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Services');
  }
}