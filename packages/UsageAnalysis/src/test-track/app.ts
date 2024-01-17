import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getIcon, getStatusIcon, Status} from './utils';
import {_package} from '../package';

interface TestCase extends Options {
  name: string;
  path: string;
  text: string;
  status: Status;
  icon: HTMLElement;
  reason: HTMLElement;
}

interface Category extends Options {
  name: string;
  children: (TestCase | Category)[];
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
  }

  async init() {
    if (this.inited) {
      grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
      return;
    }

    // Generate tree
    const filesP = _package.files.list('Test Track', true);
    const history: DG.DataFrame = await grok.functions.call('UsageAnalysis:TestTrack');
    for (const row of history.rows) {
      const path = row.get('path');
      const status: Status = row.get('status');
      const reason: string = row.get('reason');
      const el: TestCase = {name: '', path, text: '', status,
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

    // Ribbon
    const plus = ui.button(getIcon('plus', {style: 'fas'}), () => {
    }, 'Add test case');
    plus.classList.add('tt-ribbon-button');
    const report = ui.button(getIcon('tasks', {style: 'fas'}), () => {
    }, 'Generate report');
    report.classList.add('tt-ribbon-button');
    const ec = ui.button(getIcon('sort', {style: 'fas'}), () => {
      this.expanded = !this.expanded;
      this.tree.items.forEach((n: DG.TreeViewGroup | DG.TreeViewNode) => {
        if (n.constructor === DG.TreeViewGroup)
          n.expanded = this.expanded;
      });
    }, 'Expand/collapse');
    report.classList.add('tt-ribbon-button');
    const ribbon = ui.divH([plus, report, ec]);

    // Test case div
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      this.currentNode = node;
      this.testCaseDiv.innerText = node.value.text ?? '';
    });

    // UI
    this.append(ui.div([this.testCaseDiv], {id: 'tt-test-case-div-outer'}));
    this.append(ribbon);
    this.append(this.tree.root);
    this.root.style.padding = '0';
    grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
    this.inited = true;
  }

  processDir(dir: DG.FileInfo): void {
    const el: Category = {name: dir.name, children: []};
    const pathL = dir.path.split('/').slice(2);
    if (pathL.length > 1) {
      const parent = this.map[pathL.slice(0, -1).join(': ')] as Category;
      parent.children.push(el);
    } else
      this.list.push(el);
    this.map[pathL.join(': ')] = el;
  }

  async processFile(file: DG.FileInfo): Promise<void> {
    const pathL = file.path.replace(/\.[^.]+$/, '').split('/').slice(2);
    if (pathL.length < 2)
      grok.shell.error('Root test case');
    const parent = this.map[pathL.slice(0, -1).join(': ')] as Category;
    const [text, jsonS] = (await _package.files.readAsText(file)).split('\r\n---\r\n', 2);
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
      if (a.order === undefined && b.order === undefined) return 0;
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
      node.captionLabel.after(node.value.reason);
      node.captionLabel.after(node.value.icon);
      return;
    }
    const group = parent.getOrCreateGroup(obj.name, obj, false);
    for (const child of obj.children)
      this.initTreeGroupRecursive(child, group);
  }

  setContextMenu(node: DG.TreeViewNode) {
    node.captionLabel.addEventListener('contextmenu', (e) => {
      DG.Menu.popup()
        .group('Status').items(['Passed', 'Failed', 'Skipped', 'Empty'],
          (i) => {
            const status = (i === 'Empty' ? null : i.toLowerCase()) as Status;
            if (node.value.status === status) return;
            if (status === 'failed' || status === 'skipped')
              this.showChangeNodeStatusDialog(node, status);
            else
              this.changeNodeStatus(node, status);
          },
          {radioGroup: 'Status', isChecked: (i) => i.toLowerCase() === (node.value.status ?? 'empty')})
        .endGroup()
        .show();
      e.preventDefault();
      e.stopPropagation();
    });
  }

  changeNodeStatus(node: DG.TreeViewNode, status: Status, reason?: string): void {
    node.value.status = status;
    node.value.icon.innerHTML = '';
    node.value.reason.innerText = '';
    if (!status) return;
    const icon = getStatusIcon(status);
    node.value.icon.append(icon);
    if (status === 'failed' || status === 'skipped')
      node.value.reason.innerText = reason;
    const params = {success: status === 'passed', result: reason ?? '', skipped: status === 'skipped',
      type: 'manual', category: node.value.path.replace(/:\s[^:]+$/, ''), test: node.text};
    grok.log.usage(node.value.path, params, `test-manual ${node.value.path}`);
  }

  showChangeNodeStatusDialog(node: DG.TreeViewNode, status: 'failed' | 'skipped'): void {
    const dialog = ui.dialog(status === 'failed' ? 'Specify ticket' : 'Specify skip reason');
    const input = ui.textInput(status === 'failed' ? 'Key' : 'Reason', '', () => {});
    input.nullable = false;
    dialog.add(input);
    dialog.onOK(() => this.changeNodeStatus(node, status, input.value));
    dialog.show({resizable: true});
  }
};
