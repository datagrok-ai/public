import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import './../template-panel/template-panel-and-mapping.css';

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
  cancelBtn.classList.add('assay_plates__mapping-cancel-icon');
  const rightCell = ui.divH([colChoice.root, cancelBtn], 'assay_plates__dynamic-row-right-cell');
  return ui.divH([propInput.root, rightCell], 'assay_plates__mapping-editor-row');
}

export function renderMappingEditor(host: HTMLElement, options: MappingEditorOptions): void {
  const {targetProperties, sourceColumns, mappings, onMap, onUndo} = options;

  ui.empty(host);
  host.className = 'assay_plates__mapping-editor';

  if (sourceColumns.length === 0) {
    host.appendChild(ui.divText('Import a data file to map columns.', 'assay_plates__info-message'));
    return;
  }

  const tableHost = ui.divV([], 'assay_plates__mapping-editor-table');
  host.appendChild(tableHost);

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
      asterisk.className = 'assay_plates__required-asterisk';
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

    const rightCell = ui.divH([choiceControl.root], 'assay_plates__mapping-input-container');
    if (mappedSource) {
      const undoIcon = ui.iconFA('times', () => onUndo(prop.name), 'Undo mapping');
      undoIcon.classList.add('assay_plates__mapping-undo-icon');
      rightCell.appendChild(undoIcon);
    }

    const row = ui.divH([propNameEl, rightCell], 'assay_plates__mapping-editor-row');
    tableHost.appendChild(row);
  });

  const addRowHost = ui.divH([], 'assay_plates__mapping-add-row');
  const addIcon = ui.iconFA('plus', () => {
    const newDynamicRow = createDynamicMappingRow(sourceColumns, onMap, () => newDynamicRow.remove());
    tableHost.insertBefore(newDynamicRow, addRowHost);
    (newDynamicRow.querySelector('input[type="text"]') as HTMLElement)?.focus();
  }, 'Add new property mapping');
  addIcon.classList.add('assay_plates__mapping-add-icon');
  addRowHost.appendChild(addIcon);
  tableHost.appendChild(addRowHost);
}
