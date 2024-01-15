import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Octokit} from 'octokit';
import {getIcon} from './utils';
import {decode} from 'js-base64';

const octokit = new Octokit();

interface TestCase {
  name: string;
  text: string;
}

interface Category {
  name: string;
  children: (TestCase | Category)[];
}

interface OctokitResponse {
  // size: number;
  name: string;
  path: string;
  type: 'dir' | 'file' | 'submodule' | 'symlink';
  content?: string;
  // sha: string;
  // url: string;
  // git_url: string | null;
  // html_url: string | null;
  // download_url: string | null;
  // _links: {};
}

export class TestTrack extends DG.ViewBase {
  private static instance: TestTrack;
  tree: DG.TreeViewGroup;
  inited: boolean = false;
  testCaseDiv: HTMLDivElement;
  currentNode: DG.TreeViewNode | DG.TreeViewGroup;
  objTree: Category[] = [];

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
    const response = await octokit.rest.repos.getContent({
      owner: 'datagrok-ai',
      repo: 'public',
      path: 'packages/UsageAnalysis/files',
    });
    const dirs = (response.data as OctokitResponse[]).filter((f) => f.type === 'dir');
    for (const dir of dirs)
      this.objTree.push(await this.processDirRecursive(dir));
    this.objTree.forEach((obj) => this.initTreeGroupRecursive(obj, this.tree));

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

  async processDirRecursive(dir: OctokitResponse): Promise<Category> {
    const c: Category = {name: dir.name, children: []};
    const data = (await octokit.rest.repos.getContent({
      owner: 'datagrok-ai',
      repo: 'public',
      path: dir.path,
    })).data as OctokitResponse[];
    for (const child of data) {
      let res: Category | TestCase;
      if (child.type === 'dir')
        res = await this.processDirRecursive(child);
      else if (child.type === 'file')
        res = await this.processFile(child);
      else continue;
      c.children.push(res);
    }
    return c;
  }

  async processFile(file: OctokitResponse): Promise<TestCase> {
    const tc: TestCase = {name: file.name.replace(/\.[^/.]+$/, ''), text: ''};
    const data = (await octokit.rest.repos.getContent({
      owner: 'datagrok-ai',
      repo: 'public',
      path: file.path,
    })).data as OctokitResponse;
    if (data.content)
      tc.text = decode(data.content);
    return tc;
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
