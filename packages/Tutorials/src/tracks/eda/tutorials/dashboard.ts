import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {filter, tap} from 'rxjs/operators';
import {Tutorial, TutorialPrerequisites} from '@datagrok-libraries/tutorials/src/tutorial';
import {merge, of, timer} from 'rxjs';
import {elementClick} from './utils';

export class DashboardTutorial extends Tutorial {
  get name(): string {
    return 'Dashboards';
  }
  get description(): string {
    return 'Creation of interactive dashboards';
  }
  get steps(): number {
    return 26;
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
    const databasesNode = grok.shell.browsePanel.mainTree.children
      .find((child) => child.text === 'Databases') as DG.TreeViewGroup;
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
    await this.dlgInputAction(dlg, 'Set "Password" to "KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF"', 'Password',
      'KKfIh6ooS7vjzHYrNiRrderyz3KUyglrhSJF');
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
    const queryDescription = `Enter the query:<br><b>${paramAnnotation}</b><br><b>${query}</b><br>` +
      `The query uses one parameter identifying the state, and the first line, <b>${paramAnnotation}</b>, ` +
      'is its annotation. ' +
      'With it, we will be able to pass a state name into the query. ' +
      'Your query should consist of two lines: a comment with parameter ' +
      'annotation and the "select" statement that makes use of this parameter.';

    const editorText = (): string => ($(dqv.root).find('.CodeMirror')[0] as any)?.CodeMirror?.getValue() ??
      Tutorial.codeMirrorView(dqv.root)?.state.doc.toString() ?? '';
    await this.action('Enter the query with its parameter annotation',
      timer(0, 200).pipe(filter(() => {
        const lines = editorText().split('\n').map((line) => line.trim());
        return lines[0] === paramAnnotation && lines.slice(1).includes(query);
      })), () => $(dqv.root).find('.CodeMirror, .cm-editor')[0] ?? null, queryDescription,
      () => Tutorial.setCodeEditorText(dqv.root, `${paramAnnotation}\n${query}`));

    await this.buttonClickAction(
      (grok.shell.windows.simpleMode ? $('.layout-dockarea .d4-ribbon') : $('.d4-ribbon'))[0]!,
      'Save the query', 'SAVE');

    const browseSidebar = grok.shell.sidebar.getPane('Browse').header;
    await this.action('Find Browse on the sidebar and click', elementClick(() => browseSidebar), browseSidebar, '',
      () => browseSidebar.click());

    const paramEditorDlg = await this.openDialog(
      'Find the created query in the browse view, right-click it and hit Run',
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
    const uploadProjectInfo = 'Click on the "Save" button in the ribbon of the table view.';

    const projectName = 'Coffee sales dashboard';
    // every run of this tutorial leaves a dashboard with this name behind, so on a shared server
    // matching by name alone lets the last steps pick up somebody else's project
    let savedProjectId = '';
    const projectSaved = grok.events.onProjectUploaded.pipe(
      filter((p: DG.Project) => p?.friendlyName === projectName), tap((p: DG.Project) => savedProjectId = p.id));
    const isSavedProject = (p: DG.Project) => p.id === savedProjectId;

    await this.openDialog('Save a project', 'Save project', projectPaneHints, uploadProjectInfo,
      () => projectPaneHints[0].click());

    const findDialog = (title: string) => DG.Dialog.getOpenDialogs().find((dlg) => dlg.title === title);
    const nameInput = () =>
      findDialog('Save project')?.root.querySelector<HTMLInputElement>('.ui-input-editor#name') ?? null;
    const okButton = () => findDialog('Save project')?.getButton('OK') ?? null;

    await this.action(`Set the project name to "${projectName}"`,
      merge(timer(0, 200).pipe(filter(() => nameInput()?.value === projectName)), projectSaved),
      nameInput, 'Set the name before you click "OK".',
      () => nameInput() ? Tutorial.setInputValue(nameInput()!, projectName) : projectPaneHints[0].click());

    await this.action('Click "OK"', savedProjectId ? of(true) : projectSaved, okButton, '',
      () => okButton()?.click());

    const sharingDescription = 'You can share a newly created project with other users of the platform. Also, ' +
      'there is a link your project will be available at. Copy it, if you prefer this way of sharing.';
    const shareDlg = findDialog(`Share ${projectName}`);
    await this.action('Skip the sharing step', shareDlg?.onClose ?? of(true), null, sharingDescription,
      () => shareDlg?.close());

    const closeProjectDescription =
      'You can close the project by right-clicking on the sidebar and clicking "Close all"';
    await this.action('Close the project', grok.events.onProjectClosed.pipe(filter(isSavedProject)), null,
      closeProjectDescription,
      () => grok.shell.closeAll());

    // the tree is rebuilt after Close all, so the row is looked up on every poll; a click anywhere in it counts
    const dashboardsNode = () => ($('div.d4-tree-view-item-label')
      .filter((idx, el) => (el.textContent ?? '')?.startsWith('Dashboards'))[0]
      ?.closest('.d4-tree-view-node') as HTMLElement) ?? null;

    await this.action('Open browse and click on Dashboards', elementClick(dashboardsNode), dashboardsNode, '',
      () => dashboardsNode()!.click());

    await this.action('Find and open your project',
      grok.events.onProjectOpened.pipe(filter(isSavedProject)), null,
      'If the list holds more than one dashboard with this name, open the one you have just saved — ' +
      'the others belong to different runs of this tutorial and query connections you might not have access to.',
      async () => (await grok.dapi.projects.find(savedProjectId)).open());
    this.showToolbox();
    const toolbox = await Tutorial.waitFor(() => $('.d4-toolbox')[1]);

    await this.textInpAction(toolbox!, 'Set State to LA', 'State', 'LA');

    await this.buttonClickAction(toolbox!, 'Click REFRESH button', 'REFRESH');
  }
}
