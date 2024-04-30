import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {DetailedLog} from "datagrok-api/dg";

export class ReportsWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor() {
    super(ui.box());
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Top reports');
    this.order = super.addProperty('order', DG.TYPE.STRING, '2');
    const currentUserId = DG.User.current().id;
    this.root.appendChild(ui.waitBox(async () => {
      const result: DG.DataFrame = await grok.functions.call('UsageAnalysis:ReportsTop20', {'packageOwnerId': currentUserId});
      const users: {[_: string]: any} = {};
      (await grok.dapi.users.list()).forEach((user) => {
        users[user.friendlyName] = {
          'avatar': user.picture,
          'name': user.friendlyName,
          'data': user,
        };
      });
      const items = [];
      for (let i = 0; i < result.rowCount; i++) {
        // todo: add css instead of inline styles
        const currentRow = result.row(i);
        const reporter = users[currentRow.get('reporter')];
        const clock = ui.iconFA('clock', null, currentRow.get('time'));
        clock.style.marginRight = '10px';
        const portrait = ui.tooltip.bind(DG.ObjectHandler.forEntity(reporter.data)?.renderIcon(reporter.data.dart)!, () => {
          return DG.ObjectHandler.forEntity(reporter.data)?.renderTooltip(reporter.data.dart)!;
        });
        portrait.style.marginRight = '10px';
        const desc = currentRow.get('description');
        const text = ui.divText(desc === 'Auto report' ? currentRow.get('error') : desc);
        // if (currentRow.get('package_owner') === currentUserId)
        //   text.style.fontWeight = 'bold';
        const content = ui.divH([clock, portrait, text]);
        content.style.alignItems = 'center';
        const item = ui.card(content);
        item.style.overflow = 'visible';
        item.style.width = '100%';
        item.style.marginBottom = 'unset';
        item.style.border = 'none';
        item.style.marginLeft = 'auto';
        item.addEventListener('click', (e) => {
          e.preventDefault();
          DetailedLog.showReportProperties(currentRow.get('report_id'), result, i);
        });
        items.push(item);
      }
      const d = ui.divV(items);
      d.style.padding = '10px';
      d.style.lineHeight = '1';
      return d;
    }));
  }
}
