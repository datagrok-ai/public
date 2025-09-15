import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import './mapping-editor.css';

export interface TargetProperty {
  name: string;
  required?: boolean;
}

export interface MappingEditorOptions {
  targetProperties: TargetProperty[];
  sourceColumns: string[];
  mappings: Map<string, string>;
  onMap: (target: string, source: string) => void;
  onUndo: (target: string) => void;
}

function createDynamicMappingRow(
  sourceColumns: string[],
  onMap: (target: string, source: string) => void,
  onCancel: () => void
): HTMLElement {
  let propName: string | null = null;
  let sourceCol: string | null = null;

  const applyMapping = () => {
    if (propName && sourceCol)
      onMap(propName, sourceCol);
  };

  const propInput = ui.input.string('', {placeholder: 'Property name...'});
  propInput.input.addEventListener('input', () => { propName = propInput.value; });
  propInput.input.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      applyMapping();
    }
  });

  const colChoice = ui.input.choice('', {
    value: null,
    items: [null, ...sourceColumns],
    nullable: true,
    onValueChanged: (value: string | null) => {
      sourceCol = value;
      applyMapping();
    },
  });

  const cancelBtn = ui.iconFA('times', onCancel, 'Cancel');
  cancelBtn.classList.add('mapping-cancel-icon');
  const rightCell = ui.divH([colChoice.root, cancelBtn], 'dynamic-row-right-cell');
  return ui.divH([propInput.root, rightCell], 'mapping-editor-row');
}

export function renderMappingEditor(host: HTMLElement, options: MappingEditorOptions): void {
  const {targetProperties, sourceColumns, mappings, onMap, onUndo} = options;

  ui.empty(host);
  host.className = 'mapping-editor';

  if (sourceColumns.length === 0) {
    host.appendChild(ui.divText('Import a data file to map columns.', 'info-message'));
    return;
  }

  const tableHost = ui.divV([], 'mapping-editor-table');
  host.appendChild(tableHost);

  // const header = ui.divH([
  //   ui.divText('Property', {style: {fontWeight: 'bold'}}),
  //   ui.divText('Source Column', {style: {fontWeight: 'bold'}}),
  // ], 'mapping-editor-header');
  // tableHost.appendChild(header);

  const allPropsMap = new Map<string, TargetProperty>();
  targetProperties.forEach((p) => allPropsMap.set(p.name, p));
  mappings.forEach((_, targetCol) => {
    if (!allPropsMap.has(targetCol))
      allPropsMap.set(targetCol, {name: targetCol, required: false});
  });

  const allTargetProps = Array.from(allPropsMap.values()).sort((a, b) => {
    if (a.required !== b.required) return a.required ? -1 : 1;
    return a.name.localeCompare(b.name);
  });

  allTargetProps.forEach((prop) => {
    const mappedSource = mappings.get(prop.name);
    const propNameEl = ui.divH([ui.span([prop.name])]);
    if (prop.required) {
      const asterisk = ui.element('sup');
      asterisk.className = 'required-asterisk';
      asterisk.innerText = '*';
      propNameEl.appendChild(asterisk);
    }

    const choiceControl = ui.input.choice('', {
      value: mappedSource || null,
      items: [null, ...sourceColumns],
      nullable: true,
      onValueChanged: (v: string | null) => {
        if (v) onMap(prop.name, v);
        else if (mappedSource) onUndo(prop.name);
      },
    });

    const rightCell = ui.divH([choiceControl.root], 'mapping-input-container');
    if (mappedSource) {
      const undoIcon = ui.iconFA('times', () => onUndo(prop.name), 'Undo mapping');
      undoIcon.classList.add('mapping-undo-icon');
      rightCell.appendChild(undoIcon);
    }

    const row = ui.divH([propNameEl, rightCell], 'mapping-editor-row');
    tableHost.appendChild(row);
  });

  const addRowHost = ui.divH([], 'mapping-add-row');
  const addIcon = ui.iconFA('plus', () => {
    const newDynamicRow = createDynamicMappingRow(sourceColumns, onMap, () => newDynamicRow.remove());
    tableHost.insertBefore(newDynamicRow, addRowHost);
    (newDynamicRow.querySelector('input[type="text"]') as HTMLElement)?.focus();
  }, 'Add new property mapping');
  addIcon.classList.add('mapping-add-icon');
  addRowHost.appendChild(addIcon);
  tableHost.appendChild(addRowHost);
}
