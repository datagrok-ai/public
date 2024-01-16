import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getIcon} from './utils';
import {_package} from '../package';

interface TestCase {
  name: string;
  text: string;
  order?: number;
}

interface Category {
  name: string;
  children: (TestCase | Category)[];
  order?: number;
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
    const files = await _package.files.list('Test Track', true);
    for (const file of files) {
      if (file.isDirectory)
        this.processDir(file);
      else
        await this.processFile(file);
    }
    this.list.forEach((c) => this.sortCategoryRecursive(c));
    this.list.forEach((obj) => this.initTreeGroupRecursive(obj, this.tree));

    // Ribbon
    const plus = ui.button(getIcon('plus', {style: 'fas'}), () => {
    }, 'Add test case');
    plus.classList.add('tt-ribbon-button');
    const report = ui.button(getIcon('tasks', {style: 'fas'}), () => {
    }, 'Generate report');
    report.classList.add('tt-ribbon-button');
    const ribbon = ui.divH([plus, report]);

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
    const pathL = file.path.replace(/\.[^/.]+$/, '').split('/').slice(2);
    if (pathL.length < 2)
      grok.shell.error('Root test case');
    const parent = this.map[pathL.slice(0, -1).join(': ')] as Category;
    const [text, jsonS] = (await _package.files.readAsText(file)).split('\n{', 2);
    let el: TestCase = {name: file.name.replace(/\.[^/.]+$/, ''), text};
    if (jsonS) {
      const json: Options = JSON.parse('{' + jsonS);
      el = {...el, ...json};
    }
    parent.children.push(el);
    this.map[pathL.join(': ')] = el;
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
      parent.item(obj.name, obj);
      return;
    }
    const group = parent.getOrCreateGroup(obj.name, obj, false);
    for (const child of obj.children)
      this.initTreeGroupRecursive(child, group);
  }
};
