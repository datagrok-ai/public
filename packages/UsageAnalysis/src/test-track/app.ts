import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {readDataframe, writeDataframe} from './utils';

const FILENAME = 'test-cases.csv';

export class TestTrack extends DG.ViewBase {
  private static instance: TestTrack;
  tree: DG.TreeViewGroup = ui.tree();
  testCaseDiv: HTMLDivElement = ui.div();
  inited: boolean = false;
  df: DG.DataFrame = DG.DataFrame.create(0);

  public static getInstance(): TestTrack {
    if (!TestTrack.instance)
      TestTrack.instance = new TestTrack();
    return TestTrack.instance;
  }

  constructor() {
    super();
    this.name = 'Test Track';
    this.path = '/TestTrack';
    this.tree.root.id = 'tt-tree';
    this.testCaseDiv.id = 'tt-test-case-div';
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      this.testCaseDiv.innerText = node.value.test_case;
    });
  }

  async init() {
    if (this.inited) {
      grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
      return;
    }
    this.df = await readDataframe(FILENAME);
    for (const row of this.df.rows)
      this.initNodeAndParentsRecursive(row);
    this.append(this.testCaseDiv);
    this.append(this.tree.root);
    this.root.style.padding = '0';
    grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
    this.inited = true;
  }

  initNodeAndParentsRecursive(row: DG.Row): DG.TreeViewGroup | DG.TreeViewNode {
    let node: DG.TreeViewGroup | DG.TreeViewNode;
    const name = row.get('name');
    const values = {id: row.get('id'), parent_id: row.get('parent_id'), test_case: row.get('test_case')};
    let group: DG.TreeViewGroup = this.tree;
    if (values.parent_id) {
      const parentRow = this.df.row(this.findRowById(values.parent_id));
      group = this.initNodeAndParentsRecursive(parentRow) as DG.TreeViewGroup;
    }
    if (row.get('is_group'))
      node = group.getOrCreateGroup(name, values, false);
    else
      node = group.item(name, values);
    this.setContextMenu(node);
    return node;
  }

  setContextMenu(node: DG.TreeViewGroup | DG.TreeViewNode) {
    const type = node.constructor.name;
    node.captionLabel.addEventListener('contextmenu', (e) => {
      const menu = DG.Menu.popup();
      if (type === 'TreeViewGroup') {
        menu
          .item('Add group', async () => {
            this.showAddDialog(node as DG.TreeViewGroup);
          })
          .item('Add test case', async () => {
            this.showAddDialog(node as DG.TreeViewGroup, true);
          });
      }
      menu.item('Remove', async () => {
        this.showRemoveDialog(node);
      });
      menu.show();
      e.preventDefault();
      e.stopPropagation();
    });
  }

  showAddDialog(parent: DG.TreeViewGroup, testCase: boolean = false) {
    const dialog = ui.dialog(`Add new ${testCase ? 'test case' : 'group'}`);
    const nameInput = ui.textInput('Name', '', () => {});
    const testCaseInput = ui.textInput('Test case', '', () => {});
    dialog.add(nameInput);
    if (testCase)
      dialog.add(testCaseInput);
    dialog.onOK(() => this.addNode(parent, nameInput.value, testCase ? testCaseInput.value : null));
    dialog.show();
  }

  showRemoveDialog(node: DG.TreeViewGroup | DG.TreeViewNode) {
    const isGroup = node.constructor.name === 'TreeViewGroup';
    const post = `"${node.text}" ${isGroup ? 'group' : 'test case'}`;
    const dialog = ui.dialog('Remove ' + post);
    dialog.add(ui.divText('Are you sure you want to delete ' + post + '?'));
    if (isGroup) {
      dialog.add(ui.divText('The following items will also be deleted:'));
      dialog.add(ui.list((node as DG.TreeViewGroup).items.map((n) => n.text)));
    }
    dialog.onOK(async () => await this.removeNodeAndChildren(node));
    dialog.show();
  }

  async addNode(parent: DG.TreeViewGroup, name: string, testCase: string | null): Promise<void> {
    let node: DG.TreeViewGroup | DG.TreeViewNode;
    const values = {id: crypto.randomUUID(), parent_id: parent.value.id, test_case: testCase};
    if (testCase === null)
      node = parent.getOrCreateGroup(name, values, false);
    else
      node = parent.item(name, values);
    console.log(node.value);
    this.setContextMenu(node);
    this.df.rows.addNew([name, values.id, values.parent_id, testCase === null, testCase]);
    await this.updateDf();
  }

  async removeNodeAndChildren(node: DG.TreeViewGroup | DG.TreeViewNode): Promise<void> {
    this.removeRowById(node.value.id);
    this.removeChildrenRecursive([node.value.id]);
    node.remove();
    await this.updateDf();
  }

  async removeChildrenRecursive(id: string[]): Promise<void> {
    const removedIds = this.removeRowsByParentId(id);
    if (removedIds.length)
      this.removeChildrenRecursive(removedIds);
  }

  findRowById(id: string): number {
    this.df.rows.match(`id = ${id}`).select();
    const idxs = [...this.df.selection.getSelectedIndexes()];
    if (idxs.length > 1)
      grok.shell.error(id);
    return idxs[0];
  }

  // findRowsByParentId(id: string): number[] {
  //   this.df.rows.match(`parent_id = ${id}`).select();
  //   return [...this.df.selection.getSelectedIndexes()];
  // }

  removeRowById(id: string): void {
    this.df.rows.removeAt(this.findRowById(id));
  }

  removeRowsByParentId(id: string[]): string[] {
    const ind: string[] = [];
    this.df.rows.removeWhere((r) => {
      const b = id.includes(r.get('parent_id'));
      if (b) ind.push(r.get('id'));
      return b;
    });
    return ind;
  }

  async updateDf() {
    await writeDataframe(FILENAME, this.df.toCsv());
  }
}
