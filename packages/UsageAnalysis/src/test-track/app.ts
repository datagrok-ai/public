import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import { colors, FAILED, getIcon, getStatusIcon, PASSED, SKIPPED, Status } from './utils';
import { _package } from '../package';
import { Subscription } from 'rxjs';

const NEW_TESTING = 'New Testing';

interface TestCase extends Options {
  name: string;
  path: string;
  text: HTMLElement;
  status: Status | null;
  icon: HTMLDivElement;
  history: HTMLDivElement;
  reason: HTMLDivElement;
  fullReason?: string;
  datasets: string[];
  projects: string[];
  layouts: string[];
}

interface Category extends Options {
  name: string;
  children: (TestCase | Category)[];
  status: Status | null;
  icon: HTMLElement;
}

interface Options {
  order?: number;
}

interface StatusInfo {
  'User': DG.User;
  'Date': any;
  'Version': string;
  'Reason'?: HTMLElement;
  'Batch'?: string;
}

export class TestTrack extends DG.ViewBase {
  private static instance: TestTrack;
  tree: DG.TreeViewGroup;
  inited: boolean = false;
  isInitializing: boolean = false;
  testCaseDiv: HTMLDivElement;
  currentNode: DG.TreeViewNode | DG.TreeViewGroup;
  map: { [key: string]: (Category | TestCase) } = {};
  nodes: { [key: string]: (DG.TreeViewNode<any>) } = {};

  list: Category[] = [];
  expanded: boolean = false;
  version: string = grok.shell.build.client.version;
  uid: string = DG.User.current().id;
  start: string;
  nameDiv: HTMLDivElement = ui.divText('', { id: 'tt-name' });
  testingName: string;
  searchInput: DG.InputBase = ui.input.search('');

  testDescription: Map<string, string> = new Map<string, string>();
  dataSetsToOpen: string[] = [];
  projectsToOpen: string[] = [];

  pauseReportSync: HTMLButtonElement | undefined;
  runReportSync: HTMLButtonElement | undefined;
  isReporting: boolean = false;
  onRemoveSubForCollabTestingSync?: any;
  onRemoveSubForReportiingSync?: any;

  isReportSyncActive: boolean = false;
  public static getInstance(): TestTrack {
    if (!TestTrack.instance)
      TestTrack.instance = new TestTrack();
    return TestTrack.instance;
  }

  constructor(testingName?: string) {
    super();
    this.name = 'Test Track';
    this.testingName = testingName ?? '';
    this.root.classList.add('test-track');
    this.path = '/TestTrack';
    this.tree = ui.tree();
    this.tree.root.id = 'tt-tree';
    this.testCaseDiv = ui.div();
    this.testCaseDiv.id = 'tt-test-case-div';
    this.currentNode = this.tree;

    let start = localStorage.getItem('TTState');
    if (start === null) {
      start = Date.now().toString();
      localStorage.setItem('TTState', start);
      grok.log.usage(`${this.version}_${start}_${this.uid}`,
        { name: NEW_TESTING, version: this.version, uid: this.uid, start: start }, `tt-new-testing`);
    }
    this.start = start;
  }

  init() {
    if (this.isInitializing)
      return;

    if (this.inited) {
      grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
      return;
    }

    this.isInitializing = true;
    this.initTreeView().then(() => {
      this.inited = true;
      this.isInitializing = false;
    });

    this.addCollabTestingSync();
    this.addReporterSync();
  }

  public reopen() {
    if (!this.onRemoveSubForCollabTestingSync) {
      this.addCollabTestingSync()
    }
    if (!this.onRemoveSubForReportiingSync) {
      this.addReporterSync();
      if (this.pauseReportSync)
        this.pauseReportSync.style.display = 'none';
      if (this.runReportSync) 
        this.runReportSync!.style.display = 'block';
      
      this.isReporting = false;
    }
  }

  private addReporterSync() {
    this.isReportSyncActive = true;
    this.onRemoveSubForReportiingSync = grok.shell.dockManager.onClosed.subscribe((v) => {
      if (v === this.root)
        this.isReportSyncActive = false;
      if (this.onRemoveSubForReportiingSync)
        this.onRemoveSubForReportiingSync.unsubscribe();

      this.onRemoveSubForReportiingSync = undefined
    });
  }

  private addCollabTestingSync() {
    this.onRemoveSubForCollabTestingSync = grok.shell.dockManager.onClosed.subscribe((v) => {

      const updateBatchInterval = setInterval(() => {
        this.UpdateBatchData();
      }, 10000)

      if (v === this.root) {
        clearInterval(updateBatchInterval);
        if (this.onRemoveSubForCollabTestingSync)
          this.onRemoveSubForCollabTestingSync.unsubscribe();

        this.onRemoveSubForCollabTestingSync = undefined
      }
    });
  }


  private updateReportsPrefix() {
    if (this.currentNode.value && !('children' in this.currentNode.value) && this.isReporting) {
      const currentDateTime: Date = new Date();
      DG.Logger.reportPrefix = `Test Track Report \n Test Case: ${this.currentNode.value.path} \n  Test Case selection time: ${currentDateTime.toISOString()} \n`;
    }
    else
      this.resetReportsPrefix();
  }

  private resetReportsPrefix() {
    DG.Logger.reportPrefix = '';
  }

  private async initTreeView(): Promise<void> {

    let searchInvoked = false

    //searchInput 
    this.searchInput.onChanged(() => {
      // this._searchItem(); 
      if (this.searchInput.value.length > 0) {
        this.searchTreeItems(this.searchInput.value);
        searchInvoked = true;
      }
      else if (searchInvoked) {
        this.closeTreeCategories();
        searchInvoked = false;
      }
    });

    this.searchInput.input.onkeyup = (event) => {
      if (event.key === 'Escape')
        this.searchInput.fireChanged();
    };

    // Generate tree
    const filesP = _package.files.list('Test Track', true);
    const namePFromDb: string = (await grok.functions.call('UsageAnalysis:TestingName',
      { uid: this.uid, version: this.version, start: this.start })) || '"New Testing"';
    const nameP: string = this.testingName || namePFromDb.substring(1, namePFromDb.length - 1);
    const history: DG.DataFrame = await grok.functions.call('UsageAnalysis:TestTrack',
      { batchName: `"${nameP}"` });


    for (const row of history.rows) {
      const path = (row.get('test') || '').replace('Unknown: ', '');
      if (path.length === 0)
        continue;

      const status: Status = row.get('status');
      const reason: string = row.get('reason');
      const user = await grok.dapi.users.find(row.get('uid'));
      const icon = status ? ui.div(getStatusIcon(status)) : ui.div();
      let reasonTooltipValue = ui.div(this.getReason(reason ?? ''));
      this.map[path] = {
        name: '', path, text: ui.markdown(''), status, history: ui.divH([], 'tt-history'),
        icon: icon, reason: ui.div((reason.includes("\n") ? 'list' : reason), 'tt-reason'), fullReason: reason, datasets: [], projects: [], layouts: [],
      };
      const map: StatusInfo = {
        'User': user,
        'Date': row.get('date'),
        'Version': this.version,
        'Batch': nameP
      };
      if (reason)
        map['Reason'] = reasonTooltipValue;
      ui.tooltip.bind(icon, () => ui.tableFromMap(map));
    }

    const files = await filesP;
    const p: Promise<void>[] = [];
    for (const file of files) {
      if (file.isDirectory)
        this.processDir(file);
      else
        p.push(this.processFile(file));
    }
    await Promise.all(p);
    this.list.forEach((c) => this.sortCategoryRecursive(c));
    this.list.forEach((obj) => this.initTreeGroupRecursive(obj, this.tree));
    this.tree.children.forEach((c) => this.updateGroupStatusRecursiveDown(c as DG.TreeViewGroup));
    this.setContextMenu();
    for (let item of this.tree.items)
      this.testDescription.set(item.value.name, item.value.text?.textContent || '');
    // Ribbon
    const gh = ui.button(getIcon('github', { style: 'fab' }), () => {
      window.open('https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis/files/Test Track',
        '_blank')?.focus();
    }, 'Test Track folder');
    gh.classList.add('tt-ribbon-button');
    const loadBtn = ui.button(getIcon('star-of-life', { style: 'fas' }), async () => {
      for (let dataset of this.dataSetsToOpen || []) {
        try {
          const df = (await (grok.functions.eval(`OpenServerFile("${dataset}")`)))[0];
          grok.shell.addTableView(df);
        }
        catch (e) {
          grok.shell.error("could not find dataset: " + dataset);
        }
      }

      for (let project of this.projectsToOpen || []) {
        try {
          const p = await grok.dapi.projects.find(project);
          p.open();
        }
        catch (e) {
          console.error("could not find project: " + project);
        }
      }
    }, 'Open test data');
    loadBtn.classList.add('tt-ribbon-button');
    const report = ui.button(getIcon('tasks', { style: 'fas' }), () => {
      const list: { name: string, category: string, status: Status | null, reason: string }[] = [];
      Object.values(this.map).forEach((el) => {
        if (el && 'children' in el) return;
        list.push({ name: el.name, category: el.path.replace(/:\s[^:]+$/, ''), status: el.status, reason: el.fullReason || '' });
      });
      const df = DG.DataFrame.fromObjects(list)!;
      df.getCol('status').meta.colors.setCategorical(colors);
      const tv = grok.shell.addTableView(df);
      tv.grid.sort(['category', 'name']);
      tv.name = 'Report';
    }, 'Generate report');
    report.classList.add('tt-ribbon-button');
    const ec = ui.button(getIcon('sort', { style: 'fas' }), () => {
      this.expanded = !this.expanded;
      this.tree.items.forEach((n: DG.TreeViewGroup | DG.TreeViewNode) => {
        if (n.constructor === DG.TreeViewGroup)
          n.expanded = this.expanded;
      });
    }, 'Expand/collapse');
    ec.classList.add('tt-ribbon-button');
    const refresh = ui.button(getIcon('sync-alt', { style: 'fas' }), () => this.refresh(), 'Refresh');
    refresh.classList.add('tt-ribbon-button');
    ec.classList.add('tt-ribbon-button');
    const start = ui.button(getIcon('plus', { style: 'fas' }), () => this.showStartNewTestingDialog(), 'Start new testing');
    start.classList.add('tt-ribbon-button');
    this.testingName = (await nameP) ?? NEW_TESTING;
    this.nameDiv.innerText = this.testingName;
    this.nameDiv.oncontextmenu = (e) => {
      this.showEditTestingNameDialog();
      e.preventDefault();
    };
    this.pauseReportSync = ui.button(getIcon('pause', { style: 'fas' }), () => {
      this.runReportSync!.style.display = 'block';
      this.pauseReportSync!.style.display = 'none';
      this.isReporting = false;
      this.resetReportsPrefix();
    }, "Reports synchronization is running");
    this.runReportSync = ui.button(getIcon('play', { style: 'fas' }), () => {
      this.runReportSync!.style.display = 'none';
      this.pauseReportSync!.style.display = 'block';
      this.isReporting = true;
      this.updateReportsPrefix();
    }, "Reports synchronization is paused");
    this.pauseReportSync.classList.add('tt-ribbon-button');
    this.runReportSync.classList.add('tt-ribbon-button');

    this.pauseReportSync.style.display = 'none';
    this.isReporting = false;

    const ribbon = ui.divH([gh, report, ec, refresh, loadBtn, start, this.pauseReportSync, this.runReportSync, this.nameDiv]);
    ribbon.style.flexGrow = '0';

    // Test case div
    const edit = ui.button(getIcon('edit'), () => this.editTestCase(this.currentNode), 'Edit test case');
    edit.id = 'tt-edit-button';
    edit.disabled = true;
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      if (this.currentNode.constructor === DG.TreeViewNode)
        this.currentNode.value.history.style.display = 'none';
      this.currentNode = node;
      this.testCaseDiv.innerHTML = '';
      this.updateReportsPrefix();
      if (node?.value && 'children' in node.value) {
        edit.disabled = true;
        return;
      }
      this.dataSetsToOpen = node?.value?.datasets || [];
      this.projectsToOpen = node?.value?.projects || [];

      this.testCaseDiv.append(node.value.text);
      edit.disabled = false;
      node.value.history.style.display = 'flex';
      if (node.value.history.classList.contains('processed')) return;
      node.value.history.classList.add('processed');
      grok.functions.call('UsageAnalysis:LastStatuses', { path: node.value.path, batchToExclude: this.testingName }).then(async (df: DG.DataFrame) => {
        if (!df.rowCount) return;
        const first = df.row(0);
        if (first.get('version') === this.version && first.get('uid') === this.uid && first.get('start') === this.start)
          df.rows.removeAt(0);
        else if (df.rowCount === 5)
          df.rows.removeAt(4);
        const n = Math.min(df.rowCount, 5 - node.value.history.children.length);
        for (const row of df.rows) {
          if (row.idx === n) break;
          const icon = getStatusIcon(row.get('status'));
          node.value.history.append(icon);
          const user = await grok.dapi.users.find(row.get('uid'));
          const map: StatusInfo = {
            'User': user,
            'Date': row.get('date'),
            'Version': row.get('version'),
            'Batch': row.get('batchName')
          };
          const reason = row.get('reason');
          if (reason)
            map['Reason'] = this.getReason(reason);
          ui.tooltip.bind(icon, () => ui.tableFromMap(map));
        }
      });
    });

    // UI
    this.append(ui.div(this.searchInput.root));
    this.append(ui.div([this.testCaseDiv, edit], { id: 'tt-test-case-div-outer' }));
    this.append(ribbon);
    this.append(this.tree.root);
    this.root.style.padding = '0';
    grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);

    this.isInitializing = false;
  }

  private closeTreeCategories() {
    const categoriesLists = this.tree.root.getElementsByClassName('d4-tree-view-group-host');
    const categoriesTitles = this.tree.root.getElementsByClassName('d4-tree-view-tri');
    const dom = this.tree.root.getElementsByClassName('d4-tree-view-node');

    for (let i = 0; i < dom.length; i++) {
      const item = dom[i] as HTMLElement;
      if (item)
        item.classList.remove('hidden');
    }
    for (const title of categoriesTitles)
      title.classList.remove('d4-tree-view-tri-expanded');
    for (let i = 0; i < categoriesLists.length; i++) {
      let element = categoriesLists[i] as HTMLElement;
      if (!element.parentElement?.classList.contains('d4-tree-view-root'))
        element.style.display = 'none';
    }
  }

  private searchTreeItems(stringToSearch: string) {
    const dom = this.tree.root.getElementsByClassName('d4-tree-view-node');
    const regExpToSearch = new RegExp(stringToSearch.toLowerCase());
    const listToShow: HTMLElement[] = [];

    function isFitsSearchString(stringToCheck: string, description?: string): boolean {
      var result = regExpToSearch.test(stringToCheck.toLocaleLowerCase());
      if (description) {
        result = result || regExpToSearch.test(description.toLocaleLowerCase())
      }
      return result;
    }
    for (let i = 0; i < dom.length; i++) {
      const item = dom[i] as HTMLElement;
      const foundFunc = isFitsSearchString(item.textContent?.toString() || '', (this.testDescription.get(item.textContent?.toString() || '') || '').toLocaleLowerCase());
      if (foundFunc) {
        listToShow[listToShow.length] = (item);
        item.classList.remove('hidden');
      }
      else
        item.classList.add('hidden');
    }

    for (const itemToShow of listToShow) {
      let itemHierarchy: HTMLElement = itemToShow;
      while (!itemHierarchy.classList.contains('d4-tree-view-root')) {
        if (itemHierarchy.classList.contains('d4-tree-view-group')) {
          let categoryToUpdate = itemHierarchy.getElementsByClassName('d4-tree-view-tri')[0];
          if (categoryToUpdate) {
            categoryToUpdate.classList.add('d4-tree-view-tri-expanded');
            categoryToUpdate.parentElement?.classList.remove('hidden');
          }
        }

        if (itemHierarchy.style.display === 'none')
          itemHierarchy.style.display = 'block';

        if (itemHierarchy.classList.contains('hidden'))
          itemHierarchy.classList.remove('hidden');

        itemHierarchy = itemHierarchy.parentElement!;
      }

    }
  }

  processDir(dir: DG.FileInfo): void {
    const el: Category = { name: dir.name, children: [], status: null, icon: ui.div() };
    const pathL = dir.path.split('/').slice(2);
    if (pathL.length > 1) {
      const parent = this.map[pathL.slice(0, -1).join(': ') + ' C'] as Category;
      parent.children.push(el);
    } else
      this.list.push(el);
    this.map[pathL.join(': ') + ' C'] = el;
  }

  async processFile(file: DG.FileInfo): Promise<void> {
    const pathL = file.path.replace(/\.[^.]+$/, '').split('/').slice(2);
    if (pathL.length < 2)
      grok.shell.error('Root test case');
    const parent = this.map[pathL.slice(0, -1).join(': ') + ' C'] as Category
    const [textS, jsonS] = (await _package.files.readAsText(file)).split('---', 3);
    const text = ui.markdown(textS);
    const path = pathL.join(': ');
    const name = file.name.replace(/\.[^.]+$/, '');
    const elOld: TestCase | undefined = this.map[path] as TestCase;
    let el: TestCase;
    if (elOld)
      el = { ...elOld, name, text };
    else {
      el = {
        name, path, text, status: null, history: ui.divH([], 'tt-history'),
        icon: ui.div(), reason: ui.div('', 'tt-reason'), datasets: [], projects: [], layouts: [],
      };
    }
    if (jsonS)
      el = { ...el, ...JSON.parse(jsonS) };
    parent.children.push(el);
    this.map[path] = el;
  }

  sortCategoryRecursive(cat: Category): void {
    if (!cat.children.length) return;
    const cats = cat.children.filter((c) => 'children' in c) as Category[];
    cat.children.sort((a, b) => {
      if ('children' in a && !('children' in b)) return -1;
      if ('children' in b && !('children' in a)) return 1;
      if (a.order === undefined && b.order === undefined)
        return a.name.localeCompare(b.name);
      if (a.order === undefined) return 1;
      if (b.order === undefined) return -1;
      if (a.order === b.order)
        return a.name.toLocaleLowerCase().localeCompare(b.name.toLocaleLowerCase());
      return a.order > b.order ? 1 : -1;
    });
    cats.forEach((c) => this.sortCategoryRecursive(c));
  }


  initTreeGroupRecursive(obj: Category | TestCase, parent: DG.TreeViewGroup): void {
    if ('text' in obj) {
      const node = parent.item(obj.name, obj);
      node.value.reason.oncontextmenu = (e: PointerEvent) => {
        this.showNodeDialog(node, node.value.status, true);
        e.stopImmediatePropagation();
        e.preventDefault();
      };
      node.captionLabel.after(node.value.reason);
      node.captionLabel.after(node.value.history);
      node.captionLabel.after(node.value.icon);

      this.nodes[obj.path] = node;
      return;
    }

    const group = parent.getOrCreateGroup(obj.name, obj, false);
    group.captionLabel.after(group.value.icon);
    for (const child of obj.children)
      this.initTreeGroupRecursive(child, group);
  }

  setContextMenu(): void {
    this.tree.onNodeContextMenu.subscribe((data: any) => {
      const node = data?.args?.item;
      if (node?.constructor === DG.TreeViewGroup) return;
      (data.args.menu as DG.Menu)
        .group('Status').items(['Passed', 'Failed', 'Skipped'],
          (i) => {
            const status = i.toLowerCase() as Status;
            if (node.value.status === status) return;
            data.args.menu.dart.childMenuContainer?.remove();
            data.args.menu.dart.cx?.remove();
            if (status === PASSED)
              this.changeNodeStatus(node, status);
            else
              this.showNodeDialog(node, status);
          },
          { radioGroup: 'Status', isChecked: (i) => i.toLowerCase() === (node.value.status ?? 'empty') })
        .endGroup()
        .item('Edit', () => this.editTestCase(node))
        .item('Edit Reason', () => {
          this.showNodeDialog(node, data.dart.item.value.status, true);
        }, null, {
          isEnabled: () => {
            return (data.dart.item.value.status ? (data.dart.item.value.status.toLocaleLowerCase() === PASSED ? 'Status passed' : null) : 'No Status To Change');
          }
        });
    });
  }

  editTestCase(node: DG.TreeViewNode): void {
    window.open(`https://github.com/datagrok-ai/public/edit/master/packages/UsageAnalysis/files/Test Track/${node.value.path.replaceAll(': ', '/')}.md`, '_blank')?.focus();
  }

  updateGroupStatus(group: DG.TreeViewGroup): void {
    group.value.icon.innerHTML = '';
    const statuses = group.value.children.map((el: TestCase | Category) => el.status);
    if (statuses.includes(null)) {
      group.value.status = null;
      return;
    }
    const status = statuses.includes(FAILED) ? FAILED : statuses.includes(SKIPPED) ? SKIPPED : PASSED;
    group.value.status = status;
    const icon = getStatusIcon(status);
    group.value.icon.append(icon);
  }

  updateGroupStatusRecursiveUp(group: DG.TreeViewGroup): void {
    if (group === this.tree) return;
    this.updateGroupStatus(group);
    this.updateGroupStatusRecursiveUp(group.parent as DG.TreeViewGroup);
  }

  updateGroupStatusRecursiveDown(group: DG.TreeViewGroup): void {
    group.children.forEach((c) => {
      if (c.constructor === DG.TreeViewGroup)
        this.updateGroupStatusRecursiveDown(c);
    });
    this.updateGroupStatus(group);
  }

  // To Do: fix styles
  showNodeDialog(node: DG.TreeViewNode, status: typeof FAILED | typeof SKIPPED, edit: boolean = false): void {
    const name = `${edit ? 'Edit' : 'Specify'} ${status === FAILED ? 'ticket' : 'skip reason'}`;
    const dialog = ui.dialog(name);
    dialog.root.classList.add('tt-dialog', 'tt-reason-dialog');
    const value = edit ? (node.value.fullReason) || '' : '';
    const stringInput = ui.input.string(status === FAILED ? 'Key' : 'Reason', {value: value});
    stringInput.nullable = false;
    const textInput = ui.input.textArea(status === FAILED ? 'Keys' : 'Reasons', {value: value});
    textInput.nullable = false;
    let input = stringInput;
    const tabControl = ui.tabControl({
      'String': stringInput.root,
      'List': textInput.root,
    });
    if (value.includes('\n')) {
      tabControl.currentPane = tabControl.getPane('List');
      input = textInput;
    }
    dialog.root.addEventListener('keydown', (e) => {
      if (e.key == 'Enter' && tabControl.currentPane.name === 'List')
        e.stopImmediatePropagation();
    });
    tabControl.root.style.width = 'unset';
    tabControl.header.style.marginBottom = '15px';
    tabControl.onTabChanged.subscribe((tab: DG.TabPane) => input = tab.name === 'String' ? stringInput : textInput);
    dialog.add(tabControl.root);
    dialog.onOK(() => edit ? this.changeNodeReason(node, status, input.value) : this.changeNodeStatus(node, status, input.value));
    dialog.show({ resizable: true });
    dialog.initDefaultHistory();
  }

  changeNodeStatus(node: DG.TreeViewNode, status: Status, reason?: string, uid: string = this.uid, reportData: Boolean = true): void {
    const value = node.value;
    if (value.status) {
      if (value.history.children.length === 5)
        value.history.children[4].remove();
    }
    value.status = status;
    value.icon.innerHTML = '';
    value.reason.innerHTML = '';
    const icon = getStatusIcon(status);
    value.icon.append(icon);
    if (status === FAILED || status === SKIPPED) {
      if (!reason!.includes('\n'))
        value.reason.append(this.getReason(reason!));
      else
        value.reason.append(ui.label('list'));
      value.fullReason = reason;
    }
    const params = {
      success: status === PASSED, result: reason ?? '', skipped: status === SKIPPED, type: 'manual',
      category: value.path.replace(/:\s[^:]+$/, ''), name: node.text, version: this.version, uid: uid, start: this.start, ms: 0, batchName: this.testingName
    };
    if (reportData)
      grok.shell.reportTest('manual', params);
    this.updateGroupStatusRecursiveUp(node.parent as DG.TreeViewGroup);
    let reasonTooltipValue = ui.div(this.getReason(reason ?? ''));
    grok.dapi.users.find(params['uid']).then((user) => {
      const map: StatusInfo = {
        'User': user,
        'Date': dayjs(),
        'Version': params['version'],
        'Batch': this.testingName
      };

      map['Reason'] = reasonTooltipValue;
      ui.tooltip.bind(icon, () => ui.tableFromMap(map));
    });
  }

  changeNodeReason(node: DG.TreeViewNode, status: Status, reason: string, uid: string = this.uid, reportData: Boolean = true): void {
    if (status !== node.value.status) {
      grok.shell.warning('Test case status was changed');
      return;
    }
    if (reason === node.value.reason.innerText) return;
    node.value.reason.innerHTML = '';
    node.value.fullReason = reason;
    if (!reason!.includes('\n'))
      node.value.reason.append(this.getReason(reason!));
    else
      node.value.reason.append(ui.label('list'));
    const params = {
      success: status === PASSED, result: reason, skipped: status === SKIPPED, type: 'manual',
      category: node.value.path.replace(/:\s[^:]+$/, ''), name: node.text, version: this.version, uid: this.uid, start: this.start, ms: 0, batchName: this.testingName
    };

    let reasonTooltipValue = ui.div(this.getReason(reason ?? ''));
    grok.dapi.users.find(params['uid']).then((user) => {
      const map: StatusInfo = {
        'User': user,
        'Date': dayjs(),
        'Version': params['version'],
        'Batch': this.testingName
      };

      map['Reason'] = reasonTooltipValue;
      ui.tooltip.bind(node.value.icon, () => ui.tableFromMap(map));
    });
    if (reportData)
      grok.shell.reportTest('manual', params);
  }

  showStartNewTestingDialog(): void {
    const dialog = ui.dialog('Start new testing');
    dialog.root.classList.add('tt-dialog');
    const input = ui.input.string('Name', { value: NEW_TESTING });
    input.nullable = false;
    dialog.add(ui.divText('Enter name of the new testing:'));
    dialog.add(input);
    dialog.onOK(() => {
      const start = Date.now().toString();
      this.testingName = input.value;
      localStorage.setItem('TTState', start);
      grok.log.usage(`${this.version}_${start}_${this.uid}`,
        { name: this.testingName, version: this.version, uid: this.uid, start: this.start }, `tt-new-testing`);
      this.reload();
    });
    dialog.show();
  }

  showEditTestingNameDialog(): void {
    const dialog = ui.dialog('Edit testing name');
    dialog.root.classList.add('tt-dialog');
    const input = ui.input.string('Name', { value: this.testingName });
    dialog.add(input);
    dialog.onOK(() => {
      this.testingName = input.value;
      this.nameDiv.innerText = this.testingName;
      grok.log.usage(`${this.version}_${this.start}_${this.uid}`,
        { name: this.testingName, version: this.version, uid: this.uid, start: this.start }, `tt-new-testing`);
    });
    dialog.show();
  }

  reload(): void {
    const oldRoot = this.root;
    ui.setUpdateIndicator(oldRoot);
    TestTrack.instance = new TestTrack(this.testingName);
    TestTrack.getInstance().init();
    grok.shell.dockManager.close(oldRoot);
  }

  refresh(): void {
    this.UpdateBatchData();
  }

  getReason(reason: string): HTMLElement {
    if (reason.includes('\n')) {
      // const el = ui.divText(reason, 'tt-link tt-link-list');
      // el.setAttribute('data-label', 'LIST');
      // const df = DG.DataFrame.fromColumns([DG.Column.fromList("string", "key", reason.split('\n'))]);
      // const grid = df.plot.grid();
      // return grid.root;
      const res = ui.divV(reason.split('\n').map((e) => { return ui.label(e) }))
      return res;
    }
    const jira = reason.match(/GROK-\d{1,6}\b/);
    if (jira) return this.getReasonLink(reason, 'https://reddata.atlassian.net/browse/' + jira[0], jira[0]);
    const gh1 = reason.match(/#(\d{3,5})\b/);
    if (gh1) return this.getReasonLink(reason, 'https://github.com/datagrok-ai/public/issues/' + gh1[1], gh1[0]);
    const gh2 = reason.match(/https:\/\/github\.com\/datagrok-ai\/public\/issues\/(\d{3,5})/);
    if (gh2) return this.getReasonLink(reason, reason, `#${gh2[1]}`);
    const slack = reason.includes('datagrok.slack.com');
    if (slack) return this.getReasonLink(reason, reason, 'SLACK');
    const el = ui.label(reason);
    ui.tooltip.bind(el, () => reason);
    return el;
  }

  getReasonLink(reason: string, target: string, label: string): HTMLAnchorElement {
    const link = ui.link(reason, target, reason, 'tt-link');
    link.setAttribute('data-label', label);
    return link;
  }

  setContextPanelPreview(el: HTMLAnchorElement): void {
    el.onclick = (e) => {
      if (!e.ctrlKey) return;
      const ifrm = document.createElement('iframe');
      ifrm.setAttribute('src', el.href);
      ifrm.style.width = '100%';
      ifrm.style.height = '100%';
      ifrm.style.border = 'none';
      grok.shell.o = ifrm;
      e.preventDefault();
    };
  }

  private isTestCase(instance: any): instance is TestCase {
    return (instance as TestCase) !== undefined;
  }

  async UpdateBatchData() {
    if (!this.testingName)
      return;
    const history: DG.DataFrame = await grok.functions.call('UsageAnalysis:TestTrack', { batchName: `"${this.testingName}"` });

    for (const row of history.rows) {

      const status: Status = row.get('status');
      const reason: string = row.get('reason');
      const uid: string = row.get('uid');

      const path = (row.get('test') || '').replace('Unknown: ', '');
      if (path.length === 0)
        continue;

      if (this.map[path] && this.isTestCase(this.map[path])) {
        const caseValue: TestCase = this.map[path] as TestCase;

        const node = this.nodes[path];

        if (caseValue.status !== status ||
          caseValue.fullReason !== reason) {
          if (caseValue.status === row.status) {
            this.changeNodeReason(node, status, reason, uid, false);
          }
          else {
            this.changeNodeStatus(node, status, reason, uid, false);
          }
        }
      }
    }
  }
}
