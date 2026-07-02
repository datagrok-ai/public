import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subject} from 'rxjs';
import {Chem} from '../scripts-api';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {addCopyIcon} from '../utils/ui-utils';
import {MESSAGE_MALFORMED} from '../constants';
import {calculateDescriptors, getDescriptorsTree} from '../docker/api';

const _STORAGE_NAME = 'rdkit_descriptors';
const _KEY = 'selected';
let descriptors: any;

/** FuncCall editor for the `Chem:descriptorsDocker` function: edits the live funccall inputs
 * (table, molecules, selected); the platform hosts it in a dialog and runs the call on OK. */
export class DescriptorsEditor extends DG.FuncCallEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  columnInput!: DG.InputBase<DG.Column | null>;
  treeHost: HTMLElement;
  private treeControl?: DescriptorsTreeControl;
  private inputChangedSubject: Subject<any> = new Subject<any>();
  /** History applied before the async tree finished loading; consumed by {@link initTree}. */
  private pendingHistorySelection?: string[];

  constructor(private funcCall: DG.FuncCall) {
    const root = ui.divV([]);
    super(root);
    this.tableInput = ui.input.table('Table', {value: funcCall.inputs['table'] ?? grok.shell.t, nullable: false,
      onValueChanged: () => {
        funcCall.inputs['table'] = this.tableInput.value;
        this.updateColumnInput();
        this.inputChangedSubject.next();
      }, tooltipText: 'Input data frame containing molecule column'});
    root.appendChild(this.tableInput.root);
    this.treeHost = ui.div([ui.loader()]);
    root.appendChild(this.treeHost);
    this.updateColumnInput();
    // trigger changes to init the funccall
    this.tableInput.fireChanged();
    this.initTree();
  }

  private async initTree(): Promise<void> {
    const stored = await getSelected();
    const selected: string[] = this.pendingHistorySelection ?? this.funcCall.inputs['selected'] ?? stored;
    this.pendingHistorySelection = undefined;
    this.treeControl = buildDescriptorsTreeControl(selected, () => {
      this.funcCall.inputs['selected'] = this.treeControl!.getSelected();
      this.inputChangedSubject.next();
    });
    removeChildren(this.treeHost);
    this.treeHost.appendChild(this.treeControl.root);
    this.funcCall.inputs['selected'] = this.treeControl.getSelected();
    this.inputChangedSubject.next();
  }

  updateColumnInput(): void {
    const table = this.tableInput.value;
    if (table == null)
      return;
    this.columnInput?.root?.remove();
    this.columnInput = ui.input.column('Molecules', {
      table: table, nullable: false,
      filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE,
      value: table.columns.bySemType(DG.SEMTYPE.MOLECULE) ?? undefined,
      onValueChanged: () => {
        this.funcCall.inputs['molecules'] = this.columnInput.value;
        this.inputChangedSubject.next();
      }, tooltipText: 'Column with molecules to calculate descriptors for'});
    this.columnInput.root.children[0]?.classList.add('d4-chem-descriptors-molecule-column-input');
    this.root.insertBefore(this.columnInput.root, this.treeHost ?? null);
    this.columnInput.fireChanged();
  }

  get isValid(): boolean {
    return this.tableInput.value != null && this.columnInput.value != null &&
      (this.treeControl?.getSelected().length ?? 0) > 0;
  }

  inputFor(propertyName: string): DG.InputBase {
    switch (propertyName) {
    case 'table':
      return this.tableInput;
    case 'molecules':
      return this.columnInput;
    default:
      throw new Error(`Unknown property name: ${propertyName}`);
    }
  }

  getHistoryString(): string {
    return (this.treeControl?.getSelected() ?? this.funcCall.inputs['selected'] as string[] ?? []).join(',');
  }

  loadHistoryString(history: string): void {
    const selected = (history ?? '').split(',').filter((s) => s !== '');
    // the tree is built asynchronously; if history arrives first, initTree picks it up
    if (this.treeControl)
      this.treeControl.setSelected(selected); // fires onChanged -> writes funcCall.inputs['selected']
    else
      this.pendingHistorySelection = selected;
  }

  get onInputChanged(): Observable<any> {
    return this.inputChangedSubject;
  }
}

export function addDescriptorsColsToDf(table: DG.DataFrame, colsArr: DG.Column[]) {
  for (const col of colsArr) {
    if (!col)
      continue;
    const unusedName = table.columns.getUnusedName(col.name);
    col.name = unusedName;
    table.columns.add(col);
  }
}

/** Adds descriptors to table */
export async function addDescriptors(smilesCol: DG.Column, viewTable: DG.DataFrame): Promise<void> {
  openDescriptorsDialog(await getSelected(), async (selected: any) => {
    grok.userSettings.add(_STORAGE_NAME, _KEY, JSON.stringify(selected));
    selected = await getSelected();
    const pi = DG.TaskBarProgressIndicator.create('Calculating descriptors...');
    try {
      const cols = (await calculateDescriptors(smilesCol, selected));
      addDescriptorsColsToDf(viewTable, cols);
    } catch (e: any) {
      grok.shell.error(e);
    } finally {
      pi.close();
    }
  });
}

/** Calculates descriptors for single entry*/
export function getDescriptorsSingle(smiles: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  smiles = _convertMolNotation(smiles, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, rdKitModule);
  if (smiles === MESSAGE_MALFORMED)
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  const widget = new DG.Widget(ui.div());
  const result = ui.div();
  const selectButton = ui.bigButton('SELECT', async () => {
    openDescriptorsDialog(await getSelected(), (selected: any) => {
      grok.userSettings.add(_STORAGE_NAME, _KEY, JSON.stringify(selected));
      update();
    });
  });
  selectButton.classList.add('chem-descriptors-select-btn');

  const update = () => {
    removeChildren(result);
    result.appendChild(ui.loader());
    getSelected().then((selected) => {
      const smilesCol = DG.Column.string('smiles', 1).init(() => smiles);
      calculateDescriptors(smilesCol, selected)
        .then((columns: DG.Column[]) => {
          removeChildren(result);
          if (!columns.length)
            throw Error(`Descriptors haven't been calculated for ${smiles}`);
          const map: { [_: string]: any } = {};
          const colNames = columns.map((it) => it.name);
          for (const descriptor of selected) {
            const colIdx = colNames.findIndex((it) => it === descriptor);
            if (colIdx !== -1)
              map[descriptor] = columns[colIdx].get(0);
          }
          result.appendChild(ui.tableFromMap(map));
        });
    });
  };

  addCopyIcon(result, 'Descriptors');
  widget.root.appendChild(result);
  widget.root.appendChild(selectButton);

  update();

  return widget;
}

/** demo application for descriptor calculation*/
export function getDescriptorsApp(): void {
  const defaultSmiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
  let sketcherValue = defaultSmiles;

  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  let table = DG.DataFrame.create();
  table.name = 'Descriptors';
  const view = grok.shell.addTableView(table);

  const dsDiv = ui.divV([], 'grok-prop-panel');
  dsDiv.appendChild(getDescriptorsSingle(defaultSmiles).root);

  const sketcher = grok.chem.sketcher((smiles: string, _molfile: string) => {
    sketcherValue = smiles;
    removeChildren(dsDiv);
    dsDiv.appendChild(getDescriptorsSingle(smiles).root);
  }, defaultSmiles);
  const addButton = ui.bigButton('ADD', async () => {
    getSelected().then((selected) => {
      Chem.getDescriptorsPy(
        'smiles', DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'selected',
        DG.DataFrame.fromColumns([DG.Column.fromList('string', 'selected', selected)]),
      ).then((t) => {
      //grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'smiles', selected).then(t => {
        const columnNames = table.columns.names();
        if ((table.columns.length !== selected.length + 1) || selected.some((s: any) => !columnNames.includes(s))) {
          table = DG.DataFrame.create();
          table.name = 'Descriptors';
          view.dataFrame = table;
          for (const col of t.columns.toList())
            table.columns.addNew(col.name, col.type);
        }
        table.rows.addNew(t.columns.toList().map((c: any) => c.get(0)));
      });
    });
  });
  addButton.style.marginTop = '12px';
  const skDiv = ui.divV([sketcher, addButton], 'grok-prop-panel,dlg-sketcher,pure-form');

  const skNode = view.dockManager.dock(skDiv, DG.DOCK_TYPE.RIGHT, null, 'Sketcher', 0.25);
  view.dockManager.dock(dsDiv, DG.DOCK_TYPE.DOWN, skNode, 'Descriptors', 0.5);

  grok.events.onViewRemoved.subscribe((v: any) => {
    if (v.name === view.name) {
      windows.showToolbox = true;
      windows.showHelp = true;
      windows.showProperties = true;
    }
  });
}

type onOk = (selectedDescriptors: string[]) => void;

interface DescriptorsTreeControl {
  root: HTMLElement;
  getSelected(): string[];
  setSelected(names: string[]): void;
  saveHistory(): {[_: string]: string};
  loadHistory(history: any): void;
}

/** Builds the descriptors selection tree (groups with checkboxes, All/None links, counter).
 * Requires the module-level `descriptors` tree to be loaded (call `getSelected()` first). */
function buildDescriptorsTreeControl(selected: string[], onChanged?: () => void): DescriptorsTreeControl {
  const tree = ui.tree();
  tree.root.style.width = '300px';
  tree.root.style.minWidth = 'max-content';
  tree.root.style.overflowY = 'unset';

  const groups: { [_: string]: any } = {};
  const items: DG.TreeViewNode[] = [];
  const selectedDescriptors: { [_: string]: string } = {};

  const countLabel = ui.label(`${selected.length} checked`);
  countLabel.style.marginLeft = '24px';
  countLabel.style.display = 'inline-flex';

  const getSelectedNames = () => items.filter((i) => i.checked).map((i: any) => i.value['name']);

  ui.tooltip.bind(countLabel, () => {
    const names = getSelectedNames();
    if (names.length === 0)
      return 'No descriptors selected';
    const shown = names.slice(0, 30);
    const rest = names.length - shown.length;
    return ui.divV([
      ...shown.map((n) => ui.divText(n)),
      ...(rest > 0 ? [ui.divText(`...and ${rest} more`)] : []),
    ]);
  });

  function setCount(count: number) {
    countLabel.textContent = `${count} checked`;
  }

  /** Sets the group checkbox to checked / unchecked / indeterminate based on its items. */
  const updateGroupState = (group: DG.TreeViewGroup) => {
    const checkBox = group.checkBox as HTMLInputElement | null;
    if (!checkBox)
      return;
    const checkedCount = group.items.filter((i) => i.checked).length;
    const all = group.items.length > 0 && checkedCount === group.items.length;
    const partial = checkedCount > 0 && !all;
    checkBox.indeterminate = partial;
    // for the partial state the group's own checked flag is left alone: flipping it would
    // propagate to the children and destroy the very selection we are reflecting
    if (!partial && group.checked !== all)
      group.checked = all;
  };

  const updateAllGroupStates = () => {
    for (const g of Object.values(groups))
      updateGroupState(g);
  };

  const checkAll = (val: boolean) => {
    for (const g of Object.values(groups))
      g.checked = val;
    for (const i of items)
      i.checked = val;
    updateAllGroupStates();
    setCount(val ? items.length : 0);
    onChanged?.();
  };

  const selectAll = ui.label('All', {classes: 'd4-link-label', onClick: () => checkAll(true)});
  selectAll.style.marginLeft = '6px';
  selectAll.style.marginRight = '12px';
  const selectNone = ui.label('None', {classes: 'd4-link-label', onClick: () => checkAll(false)});

  const keys = Object.keys(descriptors);
  for (const groupName of keys) {
    //in case there are no descriptors in the group - do not create tree group
    if (!descriptors[groupName] || !descriptors[groupName]['descriptors'])
      continue;
    const group = tree.group(groupName, null, false);
    group.enableCheckBox();
    groups[groupName] = group;

    group.checkBox!.onchange = (_e) => {
      setCount(items.filter((i) => i.checked).length);
      if (group.checked)
        selectedDescriptors[group.text] = group.text;
      group.items.filter((i) => {
        if (i.checked)
          selectedDescriptors[i.text] = group.text;
      });
      // a user click resolves the indeterminate state: the group is now fully on or off
      (group.checkBox as HTMLInputElement).indeterminate = false;
      onChanged?.();
    };

    for (const descriptor of descriptors[groupName]['descriptors']) {
      const item = group.item(descriptor['name'], descriptor);
      item.enableCheckBox(selected.includes(descriptor['name']));
      items.push(item);

      item.checkBox!.onchange = (_e) => {
        setCount(items.filter((i) => i.checked).length);
        if (item.checked)
          selectedDescriptors[item.text] = groupName;
        updateGroupState(group);
        onChanged?.();
      };
    }
  }
  updateAllGroupStates();

  const saveInputHistory = (): any => {
    const resultHistory: { [_: string]: any } = {};
    const descriptorNames = Object.keys(selectedDescriptors);
    for (const descriptorName of descriptorNames)
      resultHistory[descriptorName] = selectedDescriptors[descriptorName];
    return resultHistory;
  };

  const setSelected = (names: string[]): void => {
    const nameSet = new Set(names);
    for (const i of items)
      i.checked = nameSet.has((i as any).value['name']);
    updateAllGroupStates();
    setCount(items.filter((i) => i.checked).length);
    onChanged?.();
  };

  // legacy dialog history: a {descriptorName: groupName} map (group entries have key === value)
  const loadInputHistory = (history: any): void => {
    setSelected(Object.keys(history));
  };

  return {
    root: ui.divV([ui.divH([selectAll, selectNone, countLabel]), tree.root]),
    getSelected: getSelectedNames,
    setSelected: setSelected,
    saveHistory: saveInputHistory,
    loadHistory: loadInputHistory,
  };
}

//description: Open descriptors selection dialog
function openDescriptorsDialog(selected: string[], onOK: onOk): void {
  const control = buildDescriptorsTreeControl(selected);
  ui.dialog('Descriptors')
    .add(control.root)
    .onOK(() => onOK(control.getSelected()))
    .show()
    .history(
      () => control.saveHistory(),
      (x) => control.loadHistory(x),
    );
}

//description: Get selected descriptors
export async function getSelected() : Promise<string[]> {
  if (!descriptors)
    descriptors = await getDescriptorsTree();
  const str = grok.userSettings.getValue(_STORAGE_NAME, _KEY);
  let selected = (str != null && str !== '') ? JSON.parse(str) : [];
  if (selected.length === 0 && descriptors['Lipinski'] && descriptors['Lipinski']['descriptors']) {
    selected = descriptors['Lipinski']['descriptors'].slice(0, 3).map((p: any) => p['name']);
    grok.userSettings.add(_STORAGE_NAME, _KEY, JSON.stringify(selected));
  }
  return selected;
}

//description: Removes all children from node
function removeChildren(node: any): void {
  while (node.firstChild)
    node.removeChild(node.firstChild);
}
