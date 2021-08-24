/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from "datagrok-api/grok";
import {divText} from "datagrok-api/ui";

export class SystemStatusWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor() {
    super(ui.panel([], 'welcome-system-widget widget-narrow'));
    let layout = ui.divV([]);

    grok.dapi.admin.getServiceInfos().then((services) => {
      // let table = ui.table(services, (item, idx) =>
      //   [`${item.key}:`, ui.span([item.status], `grok-plugin-status-${item.status.toLowerCase()}`)])
      // this.root.appendChild(table);
      this.root.appendChild(ui.divV(services.map(
        (s) => getStatus(s))));
    });

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'System');
    this.order = super.addProperty('order', DG.TYPE.STRING, '3');
  }
}
function getStatus(d:any){
  let root = ui.div();
  let icon = ui.iconFA('');
  if (d.status == "Running"){
    icon = ui.div([ui.iconFA('check')], {style:{color:'var(--green-2)', margin:'0 10px'}})
  } else if (d.status = "Failed"){
    icon = ui.div([ui.iconFA('exclamation-circle')], {style:{color:'var(--red-3)', margin:'0 10px'}})
  }
  let status  = d.status;
  let result = ui.divText(status.charAt(0).toUpperCase() + status.slice(1));
  
  result.style.color = 'var(--grey-3)';
  root.style.margin = '2px 0';

  root.append(ui.divH([
    ui.div([d.name], {style:{minWidth:'120px'}}),
    icon,
    result
  ]))
  return root;
}