import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import '../components/mapping_editor.css';

export interface TargetProperty {
  name: string;
  type?: string;
  min?: number;
  max?: number;
  semType?: string;
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
  onCancel: () => void,
): HTMLElement {
  let propName: string | null = null;
  let sourceCol: string | null = null;

  const applyMapping = () => {
    if (propName && sourceCol)
      onMap(propName, sourceCol);
  };

  const propInput = ui.input.string('', {placeholder: 'Property name...'});
  propInput.input.addEventListener('input', () => {propName = propInput.value;});
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

export function renderMappingEditor(host: HTMLElement, options: MappingEditorOptions, df: DG.DataFrame): void {
  const {targetProperties, sourceColumns, mappings, onMap, onUndo} = options;

  ui.empty(host);
  host.className = 'mapping-editor';

  if (sourceColumns.length === 0) {
    host.appendChild(ui.divText('Import a data file to map columns.', 'info-message'));
    return;
  }

  const tableHost = ui.divV([], 'mapping-editor-table');
  host.appendChild(tableHost);

  const header = ui.divH([
    ui.divText('Status', {style: {fontWeight: 'bold'}}),
    ui.divText('Property', {style: {fontWeight: 'bold'}}),
    ui.divText('Source Column', {style: {fontWeight: 'bold'}}),
  ], 'mapping-editor-header');
  tableHost.appendChild(header);

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
      asterisk.onmouseenter = (e: any) => ui.tooltip.show('The field is required', e.clientX, e.clientY);
      asterisk.onmouseleave = (e: any) => ui.tooltip.hide();
      asterisk.className = 'required-asterisk';
      asterisk.innerText = '*';
      propNameEl.appendChild(asterisk);
    }

    const choiceControl = ui.input.choice('', {
      value: mappedSource || null,
      items: [null, ...sourceColumns],
      nullable: true,
      onValueChanged: (v: string | null) => {
        if (v) {
          onMap(prop.name, v);
          const issues = validateMapping(prop, v, df);
          ui.empty(statusCell);

          if (issues.length > 0) {
            const mainIssue = issues.find((i) => i.severity === 'error') || issues[0];

            const iconName = mainIssue.severity === 'error' ? 'times-circle' : 'exclamation-triangle';
            const iconColor = mainIssue.severity === 'error' ? 'red' : 'orange';

            const icon = ui.iconFA(iconName, undefined, mainIssue.message);
            icon.style.color = iconColor;
            statusCell.appendChild(icon);
          } else {
            const icon = ui.iconFA('check-circle', undefined, 'Valid mapping');
            icon.style.color = 'green';
            statusCell.appendChild(icon);
          }

          grok.events.fireCustomEvent('mappingValidationChanged', {hasErrors: hasMappingErrors(host)});
          // if issues.length === 0 → mapping is valid (success)
        } else if (mappedSource) onUndo(prop.name);
      },
    });

    const statusCell = ui.divH([]);

    choiceControl.fireChanged();

    const row = ui.divH([statusCell, propNameEl, choiceControl.root], 'mapping-editor-row');
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

interface ValidationIssue {
  message: string;
  severity: 'warning' | 'error';
}

function validateMapping(
  prop: TargetProperty & { min?: number; max?: number },
  v: string,
  df: DG.DataFrame,
): ValidationIssue[] {
  const issues: ValidationIssue[] = [];

  const col = df.col(v);
  if (!col) {
    issues.push({ message: 'Column not found in the dataframe', severity: 'error' });
    return issues;
  }

  if (prop.type && (col.type !== prop.type)) {
    issues.push({
      message: `Type mismatch: expected ${prop.type}, got ${col.type}`,
      severity: 'error',
    });
    return issues;
  }

  if ((prop.semType && prop.semType === DG.SEMTYPE.MOLECULE) && (col.semType !== prop.semType)) {
    issues.push({
      message: `Semantic type mismatch: expected ${prop.semType}, got ${col.semType}`,
      severity: 'error',
    });
    return issues;
  }

  if ((col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT)) {
    if (prop.min != null && col.min != null && col.min < prop.min) {
      issues.push({
        message: `Values below allowed minimum (${prop.min})`,
        severity: 'error',
      });
    }

    if (prop.max != null && col.max != null && col.max > prop.max) {
      issues.push({
        message: `Values above allowed maximum (${prop.max})`,
        severity: 'error',
      });
    }
  }

  if (col.categories?.includes(''))
    issues.push({ message: 'Contains empty values', severity: 'warning' });

  return issues;
}

function hasMappingErrors(tableHost: HTMLElement): boolean {
  return Array.from(tableHost.querySelectorAll('.mapping-editor-row')).some((row) =>
    row.querySelector('.fa-times-circle') !== null,
  );
}
