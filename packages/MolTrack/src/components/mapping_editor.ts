import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import '../../css/moltrack.css';
import {requestTitleUpdate} from '../views/registration-tab';
import {MOLTRACK_MAPPING_VALIDATION_CHANGED} from '../utils/constants';

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

export function renderMappingEditor(
  host: HTMLElement,
  options: MappingEditorOptions,
  df: DG.DataFrame,
): void {
  const {targetProperties, sourceColumns, mappings, onMap, onUndo} = options;

  ui.empty(host);

  if (sourceColumns.length === 0) {
    host.appendChild(ui.divText('Import a data file to map columns.', 'info-message'));
    return;
  }

  function getAllTargetProperties(): TargetProperty[] {
    const allPropsMap = new Map<string, TargetProperty>();
    for (const p of targetProperties) allPropsMap.set(p.name, p);
    for (const [targetCol] of mappings) {
      if (!allPropsMap.has(targetCol))
        allPropsMap.set(targetCol, {name: targetCol, required: false});
    }

    return Array.from(allPropsMap.values()).sort((a, b) => {
      if (a.required !== b.required) return a.required ? -1 : 1;
      return a.name.localeCompare(b.name);
    });
  }

  function createPropCell(prop: TargetProperty): HTMLElement {
    const propNameEl = ui.divH([ui.span([prop.name])]);
    propNameEl.classList.add('moltrack-flex-row');

    if (prop.required) {
      const asterisk = ui.element('sup');
      asterisk.innerText = '*';
      asterisk.classList.add('moltrack-required-asterisk');
      asterisk.onmouseenter = (e: MouseEvent) =>
        ui.tooltip.show('The field is required', e.clientX, e.clientY);
      asterisk.onmouseleave = () => ui.tooltip.hide();
      propNameEl.appendChild(asterisk);
    }

    return propNameEl;
  }

  function createStatusCell(): HTMLElement {
    const cell = ui.div();
    cell.classList.add('moltrack-div');
    return cell;
  }

  function updateStatusCell(statusCell: HTMLElement, issues: ValidationIssue[]): void {
    ui.empty(statusCell);

    if (issues.length === 0) {
      const icon = ui.iconFA('check-circle', undefined, 'Valid mapping');
      icon.style.color = 'green';
      icon.classList.add('moltrack-div');
      statusCell.appendChild(icon);
      return;
    }

    const mainIssue = issues.find((i) => i.severity === 'error') ?? issues[0];
    const iconName = mainIssue.severity === 'error' ? 'times-circle' : 'exclamation-triangle';
    const iconColor = mainIssue.severity === 'error' ? 'red' : 'orange';
    const icon = ui.iconFA(iconName, undefined, mainIssue.message);
    icon.style.color = iconColor;
    statusCell.appendChild(icon);
  }

  function createChoiceControl(
    prop: TargetProperty,
    mappedSource: string | undefined,
    statusCell: HTMLElement,
  ): DG.InputBase<string | null> {
    const choice = ui.input.choice('', {
      value: mappedSource ?? null,
      items: [null, ...sourceColumns],
      nullable: true,
      onValueChanged: (v: string | null) => {
        if (v) {
          onMap(prop.name, v);
          const issues = validateMapping(prop, v, df);
          updateStatusCell(statusCell, issues);
        } else if (mappedSource) {
          onUndo(prop.name);
          ui.empty(statusCell);
        }

        fireValidationEvent(host);
        requestTitleUpdate();
      },
    });

    choice.root.classList.add('moltrack-input-editor');
    return choice;
  }

  const table = ui.table(
    getAllTargetProperties(),
    (prop) => {
      const mappedSource = mappings.get(prop.name);
      const statusCell = createStatusCell();
      const propCell = createPropCell(prop);
      const choiceControl = createChoiceControl(prop, mappedSource, statusCell);

      choiceControl.fireChanged();
      return [statusCell, propCell, choiceControl.root];
    },
    ['Status', 'Property', 'Source Column'],
  );

  table.classList.add('moltrack-table');
  host.appendChild(table);
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
    issues.push({message: 'Column not found in the dataframe', severity: 'error'});
    return issues;
  }

  if (prop.type) {
    const numericTypes = [DG.TYPE.INT, DG.TYPE.FLOAT];
    const bothNumeric = numericTypes.includes(prop.type as DG.TYPE) && numericTypes.includes(col.type as DG.TYPE);

    if (!bothNumeric && col.type !== prop.type) {
      issues.push({
        message: `Type mismatch: expected ${prop.type}, got ${col.type}`,
        severity: 'error',
      });
      return issues;
    }
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
    issues.push({message: 'Contains empty values', severity: 'warning'});

  return issues;
}

function fireValidationEvent(host: HTMLElement) {
  const hasErrors = host.querySelector('.moltrack-table .fa-times-circle') !== null;
  grok.events.fireCustomEvent(MOLTRACK_MAPPING_VALIDATION_CHANGED, {hasErrors: hasErrors});
}
