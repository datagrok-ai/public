import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Octokit} from 'octokit';
import {getIcon} from './utils';

const octokit = new Octokit();

export class TestTrack extends DG.ViewBase {
  private static instance: TestTrack;
  tree: DG.TreeViewGroup;
  inited: boolean = false;
  testCaseDiv: HTMLDivElement;
  currentNode: DG.TreeViewNode | DG.TreeViewGroup;

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
    const dirs = (response.data as any[]).filter((f) => f.type === 'dir');
    // console.log(dirs);

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
      // this.testCaseDiv.innerText = node.value.testCase;
    });

    // UI
    this.append(ui.div([this.testCaseDiv], {id: 'tt-test-case-div-outer'}));
    this.append(ribbon);
    this.append(this.tree.root);
    this.root.style.padding = '0';
    grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
    this.inited = true;
  }
};
