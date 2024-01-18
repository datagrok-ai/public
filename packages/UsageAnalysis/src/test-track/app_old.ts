import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {merge} from 'rxjs';
import {readDataframe, writeDataframe, getIcon, getStatusIcon, Status, FILENAME, colors} from './utils';

interface TestCaseValues {
  id: string;
  parentId: string;
  testCase: string | null;
  status: Status;
  icon: HTMLElement;
  reason: HTMLElement;
}

export class TestTrack extends DG.ViewBase {
  private static instance: TestTrack;
  tree: DG.TreeViewGroup;
  testCaseDiv: HTMLDivElement;
  inited: boolean = false;
  df: DG.DataFrame = DG.DataFrame.create(0);
  currentNode: DG.TreeViewNode | DG.TreeViewGroup;
  saveButton: HTMLButtonElement;

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
    this.tree.value = {id: null, parent_id: null, test_case: null};
    this.tree.root.id = 'tt-tree';
    this.testCaseDiv = ui.div();
    this.testCaseDiv.id = 'tt-test-case-div';
    this.currentNode = this.tree;
    this.saveButton = ui.button(getIcon('save', {style: 'fas'}), async () => {
      this.saveButton.disabled = true;
      await this.updateDf();
    }, 'Save changes');
    this.saveButton.classList.add('tt-ribbon-button');
    this.saveButton.id = 'tt-ribbon-save-button';
    this.saveButton.disabled = true;
  }

  async init() {
    if (this.inited) {
      grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
      return;
    }

    // Generate tree
    this.df = await readDataframe(FILENAME);
    const order = [...this.df.getSortedOrder(['parent_id', 'is_group', 'name'], [true, true, false])].reverse();
    let row: DG.Row;
    for (const i of order) {
      row = this.df.row(i);
      this.initNodeAndParentsRecursive(row);
    }
    merge(this.df.onDataChanged, this.df.onRowsAdded, this.df.onRowsRemoved)
      .subscribe(() => {
        this.saveButton.disabled = false;
      });

    // Ribbon
    const plus = ui.button(getIcon('plus', {style: 'fas'}), () => {
      this.showAddDialog(this.tree);
    }, 'Add root group');
    plus.classList.add('tt-ribbon-button');
    const report = ui.button(getIcon('tasks', {style: 'fas'}), () => {
      const df1 = this.df.clone(undefined, ['parent_id', 'name', 'is_group', 'status']);
      const df2 = this.df.clone(undefined, ['id', 'name']);
      df1.join(df2, ['parent_id'], ['id'], null, null, 'left', true);
      df1.rows.removeWhere((r) => r.get('is_group'));
      df1.columns.remove('is_group');
      df1.getCol('test-cases.name').name = 'category';
      df1.getCol('status').colors.setCategorical(colors);
      const tv = grok.shell.addTableView(df1);
      tv.grid.columns.setVisible(['category', 'name', 'status']);
      tv.grid.columns.setOrder(['category', 'name', 'status']);
      tv.grid.sort(['category']);
      tv.name = 'Report';
    }, 'Generate report');
    report.classList.add('tt-ribbon-button');
    const ribbon = ui.divH([plus, this.saveButton, report]);

    // Test case div
    const edit = ui.button(getIcon('edit'), () => {
      this.testCaseDiv.setAttribute('contenteditable', 'true');
      this.testCaseDiv.focus();
      this.testCaseDiv.addEventListener('focusout', () => {
        this.testCaseDiv.setAttribute('contenteditable', 'false');
        if (this.currentNode.value.testCase === this.testCaseDiv.innerText) return;
        this.currentNode.value.testCase = this.testCaseDiv.innerText;
        this.df.set('test_case', this.findRowById(this.currentNode.value.id), this.testCaseDiv.innerText);
      }, {once: true});
    }, 'Edit test case');
    edit.id = 'tt-edit-button';
    edit.disabled = true;
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      this.currentNode = node;
      this.testCaseDiv.innerText = node.value.testCase;
      if (node.constructor.name === 'TreeViewGroup')
        edit.disabled = true;
      else
        edit.disabled = false;
    });

    // UI
    this.append(ui.div([this.testCaseDiv, edit], {id: 'tt-test-case-div-outer'}));
    this.append(ribbon);
    this.append(this.tree.root);
    this.root.style.padding = '0';
    grok.shell.dockManager.dock(this.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.3);
    this.inited = true;
  }

  initNodeAndParentsRecursive(row: DG.Row): DG.TreeViewGroup | DG.TreeViewNode {
    let node: DG.TreeViewGroup | DG.TreeViewNode;
    const name = row.get('name');
    const status: Status = row.get('status') || null;
    const values: TestCaseValues = {id: row.get('id'), parentId: row.get('parent_id'),
      testCase: row.get('test_case'), status, icon: status ? ui.div(getStatusIcon(status)) : ui.div(), reason: ui.div([], 'tt-reason')};
    let group: DG.TreeViewGroup = this.tree;
    if (values.parentId) {
      const parentRow = this.df.row(this.findRowById(values.parentId));
      group = this.initNodeAndParentsRecursive(parentRow) as DG.TreeViewGroup;
    }
    if (row.get('is_group'))
      node = group.getOrCreateGroup(name, values, false);
    else
      node = group.item(name, values);
    this.setContextMenu(node);
    node.captionLabel.after(node.value.reason);
    node.captionLabel.after(node.value.icon);
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
      } else {
        menu.group('Status')
          .items(['Passed', 'Failed', 'Skipped', 'Empty'],
            (i) => {
              const status = (i === 'Empty' ? null : i.toLowerCase()) as Status;
              if (status === 'failed' || status === 'skipped')
                this.showChangeNodeStatusDialog(node, status);
              else
                this.changeNodeStatus(node, status);
            },
            {radioGroup: 'erter', isChecked: (i) => i.toLowerCase() === (node.value.status ?? 'empty')})
          .endGroup();
      }
      menu
        .item('Rename', async () => {
          this.showRenameDialog(node);
        })
        .item('Delete', async () => {
          // TO DO?: add option to delete group without children
          this.showRemoveDialog(node);
        });
      menu.show();
      e.preventDefault();
      e.stopPropagation();
    });
  }

  changeNodeStatus(node: DG.TreeViewGroup | DG.TreeViewNode, status: Status, reason?: string): void {
    if (node.value.status === status) return;
    node.value.status = status;
    node.value.icon.innerHTML = '';
    node.value.reason.innerText = '';
    this.df.set('status', this.findRowById(node.value.id), status);
    // TO DO: save category to node.value & name recursive
    const parentRow = this.df.rows.get(this.findRowById(node.value.parentId));
    if (!status) return;
    const icon = getStatusIcon(status);
    node.value.icon.append(icon);
    if (status === 'failed' || status === 'skipped')
      node.value.reason.innerText = reason;
    const params = {success: status === 'passed', result: reason ?? '', skipped: status === 'skipped',
      type: 'manual', category: parentRow.get('name'), test: node.text};
    grok.log.usage(node.value.id, params, `test-manual ${params.category}: ${params.test}`);
  }

  showAddDialog(parent: DG.TreeViewGroup, testCase: boolean = false): void {
    const dialog = ui.dialog(`Add new ${testCase ? 'test case' : 'group'}`);
    dialog.root.addEventListener('keydown', (e) => {
      if (e.key == 'Enter') {
        e.stopImmediatePropagation();
        e.stopPropagation();
      }
    });
    const nameInput = ui.textInput('Name', '', () => {});
    const testCaseInput = ui.textInput('Test case', '', () => {});
    nameInput.nullable = false;
    dialog.add(nameInput);
    if (testCase) {
      dialog.add(testCaseInput);
      testCaseInput.input.classList.add('tt-new-test-case-input');
      testCaseInput.root.style.maxWidth = 'unset';
    }
    dialog.onOK(() => this.addNode(parent, nameInput.value, testCase ? testCaseInput.value : null));
    dialog.show({resizable: true});
  }

  showRemoveDialog(node: DG.TreeViewGroup | DG.TreeViewNode): void {
    const isGroup = node.constructor.name === 'TreeViewGroup';
    const post = `"${node.text}" ${isGroup ? 'group' : 'test case'}`;
    const dialog = ui.dialog('Remove ' + post);
    dialog.add(ui.divText('Are you sure you want to remove ' + post + '?'));
    if (isGroup) {
      dialog.add(ui.divText('The following items will also be deleted:'));
      dialog.add(ui.list((node as DG.TreeViewGroup).items.map((n) => n.text)));
    }
    dialog.onOK(async () => this.removeNodeAndChildren(node));
    dialog.show();
  }

  showRenameDialog(node: DG.TreeViewGroup | DG.TreeViewNode): void {
    const dialog = ui.dialog(`Rename "${node.text}"`);
    const nameInput = ui.textInput('New name', '', () => {});
    nameInput.nullable = false;
    dialog.add(nameInput);
    dialog.onOK(() => this.renameNode(node, nameInput.value));
    dialog.show({resizable: true});
  }

  // TO DO: improve if not gh
  showChangeNodeStatusDialog(node: DG.TreeViewGroup | DG.TreeViewNode, status: 'failed' | 'skipped'): void {
    const dialog = ui.dialog(status === 'failed' ? 'Specify task' : 'Specify skip reason');
    const input = ui.textInput(status === 'failed' ? 'Key' : 'Reason', '', () => {});
    input.nullable = false;
    dialog.add(input);
    dialog.onOK(() => this.changeNodeStatus(node, status, input.value));
    dialog.show({resizable: true});
  }

  addNode(parent: DG.TreeViewGroup, name: string, testCase: string | null): void {
    let node: DG.TreeViewGroup | DG.TreeViewNode;
    const values: TestCaseValues = {id: crypto.randomUUID(), parentId: parent.value.id,
      testCase: testCase, status: 'failed', icon: ui.div(), reason: ui.div([], 'tt-reason')};
    if (testCase === null)
      node = parent.getOrCreateGroup(name, values, false);
    else
      node = parent.item(name, values);
    this.setContextMenu(node);
    node.captionLabel.after(node.value.reason);
    node.captionLabel.after(node.value.icon);
    this.df.rows.addNew([name, values.id, values.parentId, testCase === null, testCase, null]);
  }

  renameNode(node: DG.TreeViewGroup | DG.TreeViewNode, name: string): void {
    node.captionLabel.innerText = name;
    this.df.set('name', this.findRowById(node.value.id), name);
  }

  removeNodeAndChildren(node: DG.TreeViewGroup | DG.TreeViewNode): void {
    this.removeRowById(node.value.id);
    this.removeChildrenRecursive([node.value.id]);
    node.remove();
  }

  removeChildrenRecursive(id: string[]): void {
    const removedIds = this.removeRowsByParentIds(id);
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

  removeRowsByParentIds(id: string[]): string[] {
    const ind: string[] = [];
    this.df.rows.removeWhere((r) => {
      const b = id.includes(r.get('parent_id'));
      if (b) ind.push(r.get('id'));
      return b;
    });
    return ind;
  }

  async updateDf(): Promise<void> {
    await writeDataframe(FILENAME, this.df.toCsv());
  }
}
