// shared/mapping-utils.ts
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
export interface DynamicMappingRowOptions {
  sourceColumns: string[];
  onMap: (target: string, source: string) => void;
  onCancel: () => void;
  placeholder?: string;
}

export function createDynamicMappingRow(options: DynamicMappingRowOptions): HTMLElement {
  const {sourceColumns, onMap, onCancel, placeholder = 'Property name...'} = options;

  let propName: string | null = null;
  let sourceCol: string | null = null;

  const applyMapping = () => {
    if (propName && sourceCol)
      onMap(propName, sourceCol);
  };

  const propInput = ui.input.string('', {placeholder});
  propInput.onChanged.subscribe(() => {
    propName = propInput.value;
  });

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
  cancelBtn.classList.add('assay-plates--mapping-cancel-icon');

  const rightCell = ui.divH([colChoice.root, cancelBtn], 'assay-plates--mapping-input-container');
  const row = createFormRow(propInput.root, rightCell);

  (propInput.input as HTMLElement).focus();
  return row;
}

export function createFormRow(label: string | HTMLElement, input: DG.InputBase<any> | HTMLElement): HTMLElement {
  const labelEl = typeof label === 'string' ? ui.divText(label, 'ui-label') : label;
  const inputEl = (input instanceof DG.InputBase) ? input.root : input;
  inputEl.querySelector('label')?.remove();
  return ui.divH([labelEl, inputEl], 'assay-plates--form-row');
}
