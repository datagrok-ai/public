import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial, TutorialPrerequisites } from '@datagrok-libraries/tutorials/src/tutorial';
import { interval, Observable} from 'rxjs';
import { elementClick } from './utils';

export class DashboardTutorial extends Tutorial {
  get name(): string {
    return 'Dashboards';
  }
  get description(): string {
    return 'Creation of interactive dashboards';
  }
  get steps(): number {
    return 27;
  }

  get icon(): string {
    return '📊✨';
  }

  demoTable: string = '';
  helpUrl: string = '';
  prerequisites: TutorialPrerequisites = {grokConnect: true};

  protected async _run(): Promise<void> {
    grok.shell.windows.showToolbox = true;
    this.showBrowse();
    const databasesNode = grok.shell.browsePanel.mainTree.children.find((child) => child.text === 'Databases') as DG.TreeViewGroup;
    if (databasesNode) {
      databasesNode.expanded = true;
      const postgresNode = databasesNode.children.find((child) => child.text === 'Postgres') as DG.TreeViewGroup;
      if (postgresNode)
        postgresNode.expanded = true;
    }
    this.header.textContent = this.name;
    this.describe('In this tutorial, we will learn how to query data and visualize the results.');

    this.title('Access data');

    const dbViewInfo = 'In this view, you can manage connections to various data ' +
      'providers and run data queries. Each tree branch corresponds to a provider and shows ' +
      'connections to the given data source.';

    const connectionName = 'Starbucks';
    const queryName = 'Stores in @state';

    const providerRoot = $('div.d4-tree-view-group-label').filter((idx, el) =>
      (el.textContent ?? '')?.startsWith('Postgres'))[0]!;

    const dlg = await this.openDialog('Create a connection to Postgres server', 'Add new connection',
      providerRoot, `${dbViewInfo}\nOpen the context menu on the Postgres connector and click "New connection..."`,
      () => Tutorial.runContextMenu(providerRoot, 'New connection...'));

    await this.dlgInputAction(dlg, `Set "Name" to "${connectionName}"`, 'Name', connectionName);
    await this.dlgInputAction(dlg, 'Set "Server" to "db.datagrok.ai"', 'Server', 'db.datagrok.ai');
    await this.dlgInputAction(dlg, 'Set "Port" to "54324"', 'Port', '54324');
    await this.dlgInputAction(dlg, 'Set "Db" to "starbucks"', 'Db', 'starbucks');
    await this.dlgInputAction(dlg, 'Set "Login" to "datagrok"', 'Login', 'datagrok');
    await this.dlgInputAction(dlg, 'Set "Password" to "KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF"', 'Password', 'KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF');
    await this.dialogOkAction(dlg);

    const starbucksNode = () => $(providerRoot).closest('.d4-tree-view-group').find('div.d4-tree-view-group-label')
      .filter((idx, el) => el.textContent === connectionName).last()[0] ?? null;
    const dqv = await this.openViewByType(`Create a data query to the "${connectionName}" data connection`,
      'DataQueryView', starbucksNode,
      `Open the context menu on Postgres | ${connectionName} and click "New Query..."`,
      () => Tutorial.runContextMenu(starbucksNode()!, 'New Query...'));

    await this.textInpAction(dqv.root, `Set "Name" to "${queryName}"`, 'Name', queryName);

    const query = 'select * from starbucks_us where state = @state;';
    const paramAnnotation = '--input: string state';
    const queryDescription = 'As you can see, the query uses one parameter identifying the state. ' +
      'Let\'s add an annotation for it. After that, we will be able to pass a state name into the query.';
    const paramQueryDescription = 'Your query should now consist of two lines: a comment with parameter ' +
      'annotation and the "select" statement that makes use of this parameter.';

    await this.action(`Add "${query}" to the editor`,
      interval(1000).pipe(filter(() => $(dqv.root).find('pre.CodeMirror-line > span')
        .filter((idx, el) => el.textContent?.trim() === query)[0] != null)),
      null, queryDescription, () => Tutorial.setCodeEditorText(dqv.root, query),
    );

    await this.action(`Add "${paramAnnotation}" as the first line of the query`,
      interval(1000).pipe(filter(() => {
        const lines = $(dqv.root).find('pre.CodeMirror-line > span');
        return lines.length > 0 && lines[0] !== undefined && lines[0].textContent?.trim() === paramAnnotation;
      })), null, paramQueryDescription, () => Tutorial.setCodeEditorText(dqv.root, `${paramAnnotation}\n${query}`));

    await this.buttonClickAction((grok.shell.windows.simpleMode ? $('.layout-dockarea .d4-ribbon') : $('.d4-ribbon'))[0]!,
      'Save the query', 'SAVE');

    const browseSidebar = grok.shell.sidebar.getPane('Browse').header;
    await this.action('Find Browse on the sidebar and click', elementClick(() => browseSidebar), browseSidebar, '',
      () => browseSidebar.click());

    const paramEditorDlg = await this.openDialog('Find the created query in the browse view, right-click it and hit Run',
      queryName, () => $('div.d4-tree-view-item-label')
        .filter((idx, el) => (el.textContent ?? '').includes(queryName))[0] ?? null, '',
      async () => (await grok.dapi.queries.filter(`friendlyName = "${queryName}"`).order('createdOn', true).first())
        .prepare().edit());

    await this.dlgInputAction(paramEditorDlg, 'Set state to "NY"', 'State', 'NY');

    const resultRowCount = 645;
    await this.action('Click "OK" to run the query', grok.functions.onAfterRunAction.pipe(filter((call) => {
      const res = call.outputs.get('result');
      return (call.func.name === 'StoresInState' || call.func.name.includes('StoresInState_')) &&
        res instanceof DG.DataFrame && res?.rowCount === resultRowCount;
    })), $(paramEditorDlg.root).find('button.ui-btn.ui-btn-ok')[0], '', () => paramEditorDlg.getButton('OK').click());

    this.title('Create a dashboard');

    this.showToolbox();
    await this.openPlot('bar chart', (x) => x.type === DG.VIEWER.BAR_CHART);

    const projectPaneHints = [
      $('button.ui-btn').filter((idx, el) => el.textContent?.toLowerCase() === 'save')[0]!,
    ];
    const uploadProjectInfo = 'Click on the "Save" button in the scratchpad.';

    const projectName = 'Coffee sales dashboard';
    // every run of this tutorial leaves a dashboard with this name behind, so on a shared server
    // matching by name alone lets the last steps pick up somebody else's project
    let savedProjectId = '';
    const savedProjectSub = grok.events.onProjectUploaded.subscribe((p: DG.Project) => {
      if (p?.friendlyName === projectName)
        savedProjectId = p.id;
    });
    this.onClose.subscribe(() => savedProjectSub.unsubscribe());
    const isSavedProject = (p: DG.Project) => savedProjectId ? p.id === savedProjectId : p.friendlyName === projectName;

    const projectDlg = await this.openDialog('Save a project', 'Save project', projectPaneHints, uploadProjectInfo,
      () => projectPaneHints[0].click());
    await DG.delay(1000);
    const projectNameHint = $(projectDlg.root).find('.ui-input-editor#name')[0];

    await this.action(`Set the project name to "${projectName}"`, interval(2000).pipe(filter(() => projectName ===
      (<HTMLInputElement>$(projectDlg.root).find('.ui-input-editor#name')[0])?.value)),
    projectNameHint, '', () => Tutorial.setInputValue(projectNameHint as HTMLInputElement, projectName));

    // await this.action('Enable Data sync', new Observable((subscriber: any) => {
    //   $(projectDlg.root).find('.ui-input-switch').one('click', () => subscriber.next(true));
    // }), $(projectDlg.root).find('.ui-input-switch')[1]);

    const sharingDescription = 'You can share a newly created project with other users of the platform. Also, ' +
      'there is a link your project will be available at. Copy it, if you prefer this way of sharing.';
    const shareDlg = await this.openDialog('Click "OK"', `Share ${projectName}`,
      $(projectDlg.root).find('button.ui-btn.ui-btn-ok')[2], '',
      () => $(projectDlg.root).find('button.ui-btn.ui-btn-ok')[2]!.click());
    await this.action('Skip the sharing step', shareDlg.onClose, null, sharingDescription, () => shareDlg.close());

    const closeProjectDescription = 'You can close the project by right-clicking on the sidebar and clicking "Close all"';
    await this.action('Close the project', grok.events.onProjectClosed.pipe(filter(isSavedProject)), null, closeProjectDescription,
      () => grok.shell.closeAll());

    await DG.delay(1000);
    // the tree is rebuilt after Close all, so the row is looked up on every poll; a click anywhere in it counts
    const dashboardsNode = () => ($('div.d4-tree-view-item-label')
      .filter((idx, el) => (el.textContent ?? '')?.startsWith('Dashboards'))[0]?.closest('.d4-tree-view-node') as HTMLElement) ?? null;

    await this.action('Open browse and click on Dashboards', elementClick(dashboardsNode), dashboardsNode, '',
      () => dashboardsNode()!.click());

    await this.action('Find and open your project',
      grok.events.onProjectOpened.pipe(filter(isSavedProject)), null,
      'If the list holds more than one dashboard with this name, open the one you have just saved — ' +
      'the others belong to different runs of this tutorial and query connections you might not have access to.',
      async () => (await grok.dapi.projects.find(savedProjectId)).open());
    savedProjectSub.unsubscribe();
    this.showToolbox();
    await DG.delay(1000);

    await this.textInpAction($('.d4-toolbox')[1]!, 'Set State to LA', 'State', 'LA');

    await this.buttonClickAction($('.d4-toolbox')[1]!, 'Click REFRESH button', 'REFRESH');
  }
}
