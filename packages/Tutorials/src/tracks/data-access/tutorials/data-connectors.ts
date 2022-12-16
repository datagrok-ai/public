import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class DataConnectorsTutorial extends Tutorial {
  get name() {
    return 'Data Connectors';
  }
  get description() {
    return 'Direct connection to data sources and databases using the connector server';
  }
  get steps() { return 11; }
  
  demoTable: string = '';
  helpUrl: string = '';

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

    const dlg = await this.openDialog('Create a connection to PostgreSQL server',
      'Add new connection', $('.d4-tree-view-group')
        .filter((idx, el) => {
          let el = $(el).find('.d4-tree-view-group-label')[0];
          return el?.textContent === 'PostgresDart' || el?.textContent === 'PostgreSQL';
        })[0],
      'Open the context menu on the PostgreSQL connector and click "Add connection..."');

    await this.dlgInputAction(dlg, 'Set "Name" to "Starbucks"','Name', 'Starbucks');
    await this.dlgInputAction(dlg, 'Set "server" to "localhost"','Server', 'localhost');
    await this.dlgInputAction(dlg, 'Set "db" to "starbucks"', 'Db','starbucks');
    await this.dlgInputAction(dlg, 'Set "login" to "starbucks"', 'Login', 'starbucks');
    await this.dlgInputAction(dlg, 'Set "password" to "starbucks"', 'Password', 'starbucks');

    const dqv = await this.openViewByType('Create a data query to the "Starbucks" data connection',
      'DataQueryView', null,
      'Open the context menu on PostgreSQL | Starbucks and click "Add query..."');

    // UI generation delay
    await new Promise((resolve) => setTimeout(resolve, 1500));
    await this.textInpAction(dqv.root, 'Set "Name" to "Get Starbucks US"', 'Name', 'Get Starbucks US');

    await this.action('Add "select * from starbucks_us" to the editor and hit "Play"',
      grok.functions.onAfterRunAction.pipe(filter((call) => {
        const res = call.outputs.get('GetStarbucksUS');
        return call.func.name === 'GetStarbucksUS' &&
          res instanceof DG.DataFrame && res?.rowCount === 13509;
      })), null,
    );
  }
}
