import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getIcon, getStatusIcon, Status, colors, PASSED, FAILED, SKIPPED} from './utils';
import {_package} from '../package';

interface TestCase extends Options {
  name: string;
  path: string;
  text: HTMLElement;
  status: Status | null;
  icon: HTMLElement;
  reason: HTMLElement;
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

export class TestTrack extends DG.ViewBase {
  private static instance: TestTrack;
  tree: DG.TreeViewGroup;
  inited: boolean = false;
  testCaseDiv: HTMLDivElement;
  currentNode: DG.TreeViewNode | DG.TreeViewGroup;
  map: {[key: string]: (Category | TestCase)} = {};
  list: Category[] = [];
  expanded: boolean = false;
  version: string = grok.shell.build.client.version;
  userId: string = DG.User.current().id;
  start: string;

  public static getInstance(): TestTrack {
    if (!TestTrack.instance)
      TestTrack.instance = new TestTrack();
    return TestTrack.instance;
  }

  constructor() {
    super();
    this.name = 'Test Track';
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
    }
    this.start = start;
  }

  async init() {
    if (this.inited) {
      grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
      return;
    }

    // Generate tree
    const filesP = _package.files.list('Test Track', true);
    const history: DG.DataFrame = await grok.functions.call('UsageAnalysis:TestTrack',
      {version: this.version, userId: this.userId, start: this.start});
    for (const row of history.rows) {
      const path = row.get('path');
      const status: Status = row.get('status');
      const reason: string = row.get('reason');
      const el: TestCase = {name: '', path, text: ui.markdown(''), status,
        icon: status ? ui.div(getStatusIcon(status)) : ui.div(), reason: ui.div(reason, 'tt-reason')};
      this.map[path] = el;
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

    // Ribbon
    const gh = ui.button(getIcon('github', {style: 'fab'}), () => {
      window.open('https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis/files/Test Track',
        '_blank')?.focus();
    }, 'Test Track folder');
    gh.classList.add('tt-ribbon-button');
    const report = ui.button(getIcon('tasks', {style: 'fas'}), () => {
      const list: {name: string, category: string, status: Status | null, reason: string}[] = [];
      Object.values(this.map).forEach((el) => {
        if ('children' in el) return;
        list.push({name: el.name, category: el.path.replace(/:\s[^:]+$/, ''), status: el.status, reason: el.reason.innerText});
      });
      const df = DG.DataFrame.fromObjects(list)!;
      df.getCol('status').colors.setCategorical(colors);
      const tv = grok.shell.addTableView(df);
      tv.grid.sort(['category', 'name']);
      tv.name = 'Report';
    }, 'Generate report');
    report.classList.add('tt-ribbon-button');
    const ec = ui.button(getIcon('sort', {style: 'fas'}), () => {
      this.expanded = !this.expanded;
      this.tree.items.forEach((n: DG.TreeViewGroup | DG.TreeViewNode) => {
        if (n.constructor === DG.TreeViewGroup)
          n.expanded = this.expanded;
      });
    }, 'Expand/collapse');
    ec.classList.add('tt-ribbon-button');
    const refresh = ui.button(getIcon('sync-alt', {style: 'fas'}), () => this.refresh(), 'Refresh');
    refresh.classList.add('tt-ribbon-button');
    ec.classList.add('tt-ribbon-button');
    const start = ui.button(getIcon('plus', {style: 'fas'}), () => this.showStartNewTestingDialog(), 'Start new testing');
    start.classList.add('tt-ribbon-button');
    const ribbon = ui.divH([gh, report, ec, refresh, start]);

    // Test case div
    const edit = ui.button(getIcon('edit'), () => this.editTestCase(this.currentNode), 'Edit test case');
    edit.id = 'tt-edit-button';
    edit.disabled = true;
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      this.currentNode = node;
      this.testCaseDiv.innerHTML = '';
      if (node.value.text) {
        this.testCaseDiv.append(node.value.text);
        edit.disabled = false;
      } else
        edit.disabled = true;
    });

    // UI
    this.append(ui.div([this.testCaseDiv, edit], {id: 'tt-test-case-div-outer'}));
    this.append(ribbon);
    this.append(this.tree.root);
    this.root.style.padding = '0';
    grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
    this.inited = true;
  }

  processDir(dir: DG.FileInfo): void {
    const el: Category = {name: dir.name, children: [], status: null, icon: ui.div()};
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
    const parent = this.map[pathL.slice(0, -1).join(': ') + ' C'] as Category;
    const [textS, jsonS] = (await _package.files.readAsText(file)).split('\r\n---\r\n', 2);
    const text = ui.markdown(textS);
    const path = pathL.join(': ');
    const elOld: TestCase | undefined = this.map[path] as TestCase;
    const status = elOld ? elOld.status : null;
    const reason = elOld ? elOld.reason : ui.div('', 'tt-reason');
    const icon = elOld ? elOld.icon : ui.div();
    let el: TestCase = {name: file.name.replace(/\.[^.]+$/, ''), path, text, status, reason, icon};
    if (jsonS) {
      const json: Options = JSON.parse(jsonS);
      el = {...json, ...el};
    }
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
      return a.order > b.order ? 1 : -1;
    });
    cats.forEach((c) => this.sortCategoryRecursive(c));
  }

  initTreeGroupRecursive(obj: Category | TestCase, parent: DG.TreeViewGroup): void {
    if ('text' in obj) {
      const node = parent.item(obj.name, obj);
      this.setContextMenu(node);
      (node.value.reason as HTMLElement).ondblclick = () => this.showChangeReasonDialog(node, node.value.status);
      node.captionLabel.after(node.value.reason);
      node.captionLabel.after(node.value.icon);
      return;
    }
    const group = parent.getOrCreateGroup(obj.name, obj, false);
    group.captionLabel.after(group.value.icon);
    for (const child of obj.children)
      this.initTreeGroupRecursive(child, group);
  }

  setContextMenu(node: DG.TreeViewNode): void {
    node.captionLabel.addEventListener('contextmenu', (e) => {
      DG.Menu.popup()
        .group('Status').items(['Passed', 'Failed', 'Skipped'],
          (i) => {
            const status = i.toLowerCase() as Status;
            if (node.value.status === status) return;
            if (status === FAILED || status === SKIPPED)
              this.showChangeNodeStatusDialog(node, status);
            else
              this.changeNodeStatus(node, status);
          },
          {radioGroup: 'Status', isChecked: (i) => i.toLowerCase() === (node.value.status ?? 'empty')})
        .endGroup()
        .item('Edit', () => this.editTestCase(node))
        .show();
      e.preventDefault();
      e.stopPropagation();
    });
  }

  editTestCase(node: DG.TreeViewNode): void {
    window.open(`https://github.com/datagrok-ai/public/edit/master/packages/UsageAnalysis/files/Test Track/${
      node.value.path.replaceAll(': ', '/')}.md`, '_blank')?.focus();
  }

  changeNodeStatus(node: DG.TreeViewNode, status: Status, reason?: string): void {
    node.value.status = status;
    node.value.icon.innerHTML = '';
    node.value.reason.innerText = '';
    const icon = getStatusIcon(status);
    node.value.icon.append(icon);
    if (status === FAILED || status === SKIPPED)
      node.value.reason.innerText = reason;
    const params = {success: status === PASSED, result: reason ?? '', skipped: status === SKIPPED, type: 'manual',
      category: node.value.path.replace(/:\s[^:]+$/, ''), test: node.text, version: this.version, userId: this.userId, start: this.start};
    grok.log.usage(node.value.path, params, `test-manual ${node.value.path}`);
    this.updateGroupStatusRecursiveUp(node.parent as DG.TreeViewGroup);
  }

  updateGroupStatus(group: DG.TreeViewGroup) {
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

  updateGroupStatusRecursiveUp(group: DG.TreeViewGroup) {
    if (group === this.tree) return;
    this.updateGroupStatus(group);
    this.updateGroupStatusRecursiveUp(group.parent as DG.TreeViewGroup);
  }

  updateGroupStatusRecursiveDown(group: DG.TreeViewGroup) {
    group.children.forEach((c) => {
      if (c.constructor === DG.TreeViewGroup)
        this.updateGroupStatusRecursiveDown(c);
    });
    this.updateGroupStatus(group);
  }

  showChangeNodeStatusDialog(node: DG.TreeViewNode, status: typeof FAILED | typeof SKIPPED): void {
    const dialog = ui.dialog(status === FAILED ? 'Specify ticket' : 'Specify skip reason');
    const input = ui.textInput(status === FAILED ? 'Key' : 'Reason', '', () => {});
    input.nullable = false;
    dialog.add(input);
    dialog.onOK(() => this.changeNodeStatus(node, status, input.value));
    dialog.show({resizable: true});
  }

  changeNodeReason(node: DG.TreeViewNode, reason: string, status: Status): void {
    if (status !== node.value.status) {
      grok.shell.warning('Test case status was changed');
      return;
    }
    node.value.reason.innerText = reason;
    const params = {success: status === PASSED, result: reason, skipped: status === SKIPPED, type: 'manual',
      category: node.value.path.replace(/:\s[^:]+$/, ''), test: node.text, version: this.version, userId: this.userId, start: this.start};
    grok.log.usage(node.value.path, params, `test-manual ${node.value.path}`);
  }

  showChangeReasonDialog(node: DG.TreeViewNode, status: Status): void {
    const dialog = ui.dialog(status === FAILED ? 'Edit ticket' : 'Edit skip reason');
    const input = ui.textInput(status === FAILED ? 'Key' : 'Reason', node.value.reason.innerText, () => {});
    input.nullable = false;
    dialog.add(input);
    dialog.onOK(() => this.changeNodeReason(node, input.value, status));
    dialog.show({resizable: true});
  }

  showStartNewTestingDialog(): void {
    const dialog = ui.dialog('Confirm');
    dialog.add(ui.divText('Are you sure you want to start new testing?'));
    dialog.onOK(() => {
      localStorage.setItem('TTState', Date.now().toString());
      this.refresh();
    });
    dialog.show();
  }

  refresh() {
    const oldRoot = this.root;
    ui.setUpdateIndicator(oldRoot);
    TestTrack.instance = new TestTrack();
    TestTrack.getInstance().init().then(() => grok.shell.dockManager.close(oldRoot));
  }
};
