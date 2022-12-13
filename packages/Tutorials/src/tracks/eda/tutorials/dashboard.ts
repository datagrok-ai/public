import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class DashboardTutorial extends Tutorial {
  get name(): string {
    return 'Dashboards';
  }
  get description(): string {
    return 'Creation of interactive dashboards';
  }
  get steps(): number {
    return 7;
  }

  demoTable: string = '';
  helpUrl: string = '';

  protected async _run(): Promise<void> {
    this.header.textContent = this.name;
    this.describe('');

    const dbViewInfo = 'In this view, you can manage connections to various data ' +
      'providers and run data queries. Each tree branch corresponds to a provider and shows ' +
      'connections to the given data source.';

    await this.openViewByType(
      'Find "Data | Databases" in the sidebar to open the tree of connections',
      DG.View.DATABASES,
      this.getSidebarHints('Data', DG.View.DATABASES),
      dbViewInfo
    );

    const connectionName = 'Coffee_company';
    const queryName = 'Stores in @state';
    const queryHints = [
      $('div.d4-tree-view-group-label').filter((idx, el) => el.textContent == 'PostgresDart' || el.textContent == 'PostgreSQL')[0]!,
      // $('div.d4-tree-view-group-label').filter((idx, el) => el.textContent == connectionName)[0]!,
      // $('div.d4-tree-view-group-label').filter((idx, el) => el.textContent == queryName)[0]!,
    ];

    await this.openViewByType(`Run the "${queryName}" query`,
      DG.VIEW_TYPE.TABLE_VIEW, queryHints,
      'We will use the existing query for the PostgreSQL database. To execute it, double-click ' +
      'the query label. You can also choose "Run" in the context menu or in the property panel.');

    await this.openPlot('bar chart', (x) => x.type === DG.VIEWER.BAR_CHART);

    const projectPane = grok.shell.sidebar.getPane('Projects');
    const projectPaneHints = [
      projectPane.header,
      //$('button.ui-btn').filter((idx, el) => el.firstChild?.textContent === 'Upload')[0]!,
    ];

    await this.openDialog('Save a project', 'Upload project', projectPaneHints);

    await this.openViewByType('Open the project gallery', DG.View.PROJECTS,
      this.getSidebarHints('Data', DG.View.PROJECTS));

    await this.action('Find your project',
      grok.events.onAccordionConstructed.pipe(filter((acc) =>
        acc.context instanceof DG.Project && acc.context?.friendlyName == 'dashboard')));
  }
}
