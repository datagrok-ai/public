import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { Observable } from 'rxjs';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../tutorial';


export class DataConnectorsTutorial extends Tutorial {
  get name() { return 'Data Connectors'; }
  get description() {
    return 'Direct connection to data sources and databases using the connector server';
  }
  demoTable: string = '';

  protected async _run() {
    this.title('Connect to Data');

    await this.action('Find "Data | Databases" in the sidebar to open the tree of connections.',
      grok.events.onViewAdded.pipe(filter((v: DG.View) => v.type == DG.View.DATABASES)));

    this.describe('In this view, you can create queries to data connectors from the list. Each ' +
      'tree branch corresponds to a provider and shows connections to the given data source');

    this.title('Create a new data connection');

    let dlg: DG.Dialog;

    await this.action('Create a connection to PostgreSQL server: Open the context menu on ' +
    'the PostgreSQL connector and click "Add connection..."',
      grok.events.onDialogShown.pipe(filter((modal: DG.Dialog) => {
        if (modal.title == 'Add new connection') {
          dlg = modal;
          return true;
        }
        return false;
      })));

    await this.dlgInputAction(dlg!, 'Set "Name" to "Starbucks"', 'Name', 'Starbucks');
    await this.dlgInputAction(dlg!, 'Set "server" to "localhost"', 'Server', 'localhost');
    await this.dlgInputAction(dlg!, 'Set "db" to "starbucks"', 'Db', 'starbucks');
    await this.dlgInputAction(dlg!, 'Set "login" to "starbucks"', 'Login', 'starbucks');
    await this.dlgInputAction(dlg!, 'Set "password" to "starbucks"', 'Password', 'starbucks');

    this.title('Create a new data query');

    let dqv: DG.View;
    await this.action('Create a data query to the "Starbucks" data connection: ' +
      'Open context menu on PostgreSQL | Starbucks and click "Add query..."',
      grok.events.onViewAdded.pipe(filter((v: DG.View) => {
        if (v.type === 'DataQueryView') {
          dqv = v;
          return true;
        }
        return false;
      })),
    );

    // UI generation delay
    await new Promise((resolve) => setTimeout(resolve, 1500));
    await this.textInpAction(dqv!.root, 'Set "Name" to "Get Starbucks US"', 'Name', 'Get Starbucks US');
    this.describe('Write the following text to the editor: "select * from starbucks_us"');

    await this.action('Click on the "Play" button to run this query',
      grok.functions.onBeforeRunAction.pipe(filter((call: DG.FuncCall) => call.func.name === 'GetStarbucksUS')),
    );
  }
}
