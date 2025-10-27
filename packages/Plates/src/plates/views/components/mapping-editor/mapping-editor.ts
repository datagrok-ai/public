import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {createDynamicMappingRow} from '../../shared/mapping-utils';

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

export function renderMappingEditor(host: HTMLElement, options: MappingEditorOptions): void {
  const {targetProperties, sourceColumns, mappings, onMap, onUndo} = options;

  ui.empty(host);
  host.className = 'assay-plates--mapping-editor';

  if (sourceColumns.length === 0) {
    host.appendChild(ui.divText('Import a data file to map columns.', 'assay-plates--info-message'));
    return;
  }

  function getAllTargetProperties(): TargetProperty[] {
    const allPropsMap = new Map<string, TargetProperty>();
    targetProperties.forEach((p) => allPropsMap.set(p.name, p));
    mappings.forEach((_, targetCol) => {
      if (!allPropsMap.has(targetCol)) allPropsMap.set(targetCol, {name: targetCol, required: false});
    });
    return Array.from(allPropsMap.values()).sort((a, b) => {
      if (a.required !== b.required) return a.required ? -1 : 1;
      return a.name.localeCompare(b.name);
    });
  }

  function createPropertyCell(prop: TargetProperty): HTMLElement {
    const propNameEl = ui.divH([ui.span([prop.name])]);
    if (prop.required) {
      const asterisk = ui.element('sup');
      asterisk.className = 'assay-plates--required-asterisk';
      asterisk.innerText = '*';
      propNameEl.appendChild(asterisk);
    }
    return propNameEl;
  }

  function createChoiceCell(prop: TargetProperty): HTMLElement {
    const mappedSource = mappings.get(prop.name);
    const choiceControl = ui.input.choice('', {
      value: mappedSource ?? null,
      items: [null, ...sourceColumns],
      nullable: true,
      onValueChanged: (v: string | null) => {
        if (v) onMap(prop.name, v);
        else if (mappedSource) onUndo(prop.name);
      },
    });

    const rightCell = ui.divH([choiceControl.root], 'assay-plates--mapping-input-container');

    if (mappedSource) {
      const undoIcon = ui.iconFA('times', () => onUndo(prop.name), 'Undo mapping');
      undoIcon.classList.add('assay-plates--mapping-undo-icon');
      rightCell.appendChild(undoIcon);
    }

    return rightCell;
  }

  // function createAddRowButton(): HTMLElement {
  //   const addRowHost = ui.divH([], 'assay-plates--mapping-add-row');
  //   const addIcon = ui.iconFA('plus', () => {
  //     const newDynamicRow = createDynamicMappingRow({
  //       sourceColumns,
  //       onMap,
  //       onCancel: () => newDynamicRow.remove(),
  //     });
  //     table.parentElement?.insertBefore(newDynamicRow, addRowHost);
  //     (newDynamicRow.querySelector('input[type="text"]') as HTMLElement)?.focus();
  //   }, 'Add new property mapping');

  //   addIcon.classList.add('assay-plates--mapping-add-icon');
  //   addRowHost.appendChild(addIcon);
  //   return addRowHost;
  // }

  const table = ui.table(
    getAllTargetProperties(),
    (prop) => [createPropertyCell(prop), createChoiceCell(prop)],
    ['Property', 'Source Column']
  );

  table.classList.add('assay-plates--mapping-editor-table');
  host.appendChild(table);
  // host.appendChild(createAddRowButton());
}
