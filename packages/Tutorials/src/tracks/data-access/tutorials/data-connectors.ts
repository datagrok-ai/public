import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial, TutorialPrerequisites } from '@datagrok-libraries/tutorials/src/tutorial';


export class DataConnectorsTutorial extends Tutorial {
  get name() {
    return 'Data Connectors';
  }
  get description() {
    return 'Direct connection to data sources and databases using the connector server';
  }
  get steps() { return 13; }
  
  demoTable: string = '';
  helpUrl: string = '';
  prerequisites: TutorialPrerequisites = {grokConnect: true};

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('In this tutorial, we will browse the tree of connections and learn how to ' +
      'create new connections to query the database.');
    
    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    const dbViewInfo = 'In this view, you can create queries to data connectors from the list. ' +
      'Each tree branch corresponds to a provider and shows connections to the given data source.';

    await this.openViewByType(
      'Find "Data | Databases" in the sidebar to open the tree of connections',
      DG.View.DATABASES, this.getSidebarHints('Data', DG.View.DATABASES), dbViewInfo
    );

    const providerRoot = $('.d4-tree-view-group').filter((idx, el) => {
      const label = $(el).find('.d4-tree-view-group-label')[0];
      return (label?.textContent ?? '').startsWith('Postgres');
    })[0];

    const dlg = await this.openDialog('Create a connection to PostgreSQL server', 'Add new connection',
      providerRoot, 'Open the context menu on the PostgreSQL connector and click "Add connection..."');

    await this.dlgInputAction(dlg, 'Set "Name" to "Starbucks"', 'Name', 'Starbucks');
    await this.dlgInputAction(dlg, 'Set "Server" to "db.datagrok.ai"', 'Server', 'db.datagrok.ai');
    await this.dlgInputAction(dlg, 'Set "Port" to "54324"', 'Port', '54324');
    await this.dlgInputAction(dlg, 'Set "Db" to "starbucks"', 'Db','starbucks');
    await this.dlgInputAction(dlg, 'Set "Login" to "datagrok"', 'Login', 'datagrok');
    await this.dlgInputAction(dlg, 'Set "Password" to "datagrok"', 'Password', 'datagrok');
    await this.action('Click "OK"', dlg.onClose, $(dlg.root).find('button.ui-btn.ui-btn-ok')[0]);

    const dqv = await this.openViewByType('Create a data query to the "Starbucks" data connection',
      'DataQueryView', $(providerRoot).find('div.d4-tree-view-group-label').filter((idx, el) =>
        el.textContent === 'Starbucks')[0],
      'Open the context menu on PostgreSQL | Starbucks and click "Add query..."');

    // UI generation delay
    await new Promise((resolve) => setTimeout(resolve, 1500));
    await this.textInpAction(dqv.root, 'Set "Name" to "Get Starbucks US"', 'Name', 'Get Starbucks US');

    await this.action('Add "select * from starbucks_us" to the editor and hit "Play"',
      grok.functions.onAfterRunAction.pipe(filter((call) => {
        const res = call.outputs.get('result');
        return call.func.name === 'GetStarbucksUS' &&
          res instanceof DG.DataFrame && res?.rowCount === 13509;
      })), null,
    );
  }
}
