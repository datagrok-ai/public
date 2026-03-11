import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {KnimeInputParam, KnimeParamType, KnimeExecutionInput} from './types';
import {dataFrameToKnimeTable} from './data-conversion';

export interface WorkflowFormResult {
  root: HTMLElement;
  getInputs: () => Promise<KnimeExecutionInput>;
}

const knimeTypeToDgInputType: Record<KnimeParamType, string> = {
  'table': DG.InputType.Table,
  'string': DG.InputType.Text,
  'int': DG.InputType.Int,
  'double': DG.InputType.Float,
  'boolean': DG.InputType.Bool,
  'file': DG.InputType.File,
  'json': DG.InputType.TextArea,
};

/** Build a dynamic input form based on KNIME workflow parameter descriptors. */
export function buildWorkflowInputForm(params: KnimeInputParam[]): WorkflowFormResult {
  const inputs: {param: KnimeInputParam; input: DG.InputBase}[] = [];

  // Separate ungrouped (top-level) and grouped params
  const ungrouped: KnimeInputParam[] = [];
  const grouped = new Map<string, KnimeInputParam[]>();
  for (const param of params) {
    if (param.group) {
      if (!grouped.has(param.group))
        grouped.set(param.group, []);
      grouped.get(param.group)!.push(param);
    }
    else
      ungrouped.push(param);
  }

  const formElements: HTMLElement[] = [];

  const createInfoIcon = (tooltipContent: HTMLElement) => {
    const infoIcon = ui.iconFA('info-circle', null);
    infoIcon.classList.add('knime-info-icon');
    ui.tooltip.bind(infoIcon, () => tooltipContent);
    return infoIcon;
  }

  // Ungrouped inputs first
  if (ungrouped.length > 0) {
    const ungroupedInputs: DG.InputBase[] = [];
    for (const param of ungrouped) {
      const input = createInput(param);
      inputs.push({param, input});
      ungroupedInputs.push(input);

      const tooltipContent = buildTooltipContent(param);
      if (tooltipContent) {
        const infoIcon = createInfoIcon(tooltipContent);
        input.root.append(infoIcon);
      }
    }
    formElements.push(ui.inputs(ungroupedInputs));
  }

  // Grouped inputs as accordion panels
  if (grouped.size > 0) {
    const acc = ui.accordion();
    for (const [groupName, groupParams] of grouped) {
      const pane = acc.addPane(groupName, () => {
        const groupInputs: DG.InputBase[] = [];
        for (const param of groupParams) {
          const input = createInput(param);
          inputs.push({param, input});
          groupInputs.push(input);
        }
        return ui.inputs(groupInputs);
      }, true);

      const groupDesc = groupParams[0]?.groupDescription;
      if (groupDesc) {
        const infoIcon = createInfoIcon(ui.divText(groupDesc));
        pane.root.querySelector('.d4-accordion-pane-header')?.appendChild(infoIcon);
      }
    }
    formElements.push(acc.root);
  }

  const form = ui.div(formElements, 'knime-input-form');

  return {
    root: form,
    getInputs: async () => {
      const result: KnimeExecutionInput = {};

      for (const {param, input} of inputs) {
        const val = input.value;
        if (val === null || val === undefined)
          continue;

        let converted: any;
        if (param.type === 'table')
          converted = dataFrameToKnimeTable(val as DG.DataFrame);
        else if (param.type === 'file') {
          const fileInfo = val as DG.FileInfo;
          const bytes = await fileInfo.readAsBytes();
          converted = new File([bytes.buffer as ArrayBuffer], fileInfo.name);
        }
        else if (param.type === 'json')
          converted = val ? JSON.parse(val as string) : null;
        else
          converted = val;

        if (param.group) {
          if (!result[param.group])
            result[param.group] = {};
          result[param.group][param.name] = converted;
        }
        else
          result[param.name] = converted;
      }
      return result;
    },
  };
}

function createInput(param: KnimeInputParam): DG.InputBase {
  const dgType = knimeTypeToDgInputType[param.type] ?? DG.InputType.Text;
  const input = ui.input.forInputType(dgType);
  input.caption = param.name;
  input.nullable = !param.required;

  if (param.type === 'table') {
    (input as any).items = grok.shell.tables;
    if (grok.shell.tables.length > 0)
      input.value = grok.shell.tables[0];
  }
  else if (param.type === 'json')
    input.value = param.defaultValue ? JSON.stringify(param.defaultValue, null, 2) : '';
  else if (param.defaultValue !== undefined)
    input.value = param.defaultValue;

  return input;
}

function buildTooltipContent(param: KnimeInputParam): HTMLElement | null {
  const parts: HTMLElement[] = [];

  if (param.description)
    parts.push(ui.divText(param.description));

  if (param.tableSpec && param.tableSpec.length > 0) {
    const rows = param.tableSpec.map((col) =>
      `<tr><td>${col.name}</td><td>${col.type}</td></tr>`,
    ).join('');
    const tableHtml = ui.div([]);
    tableHtml.innerHTML =
      `<div style="margin-top:4px"><b>Table schema:</b></div>` +
      `<table class="knime-schema-table"><tr><th>Column</th><th>Type</th></tr>${rows}</table>`;
    parts.push(tableHtml);
  }

  if (parts.length === 0)
    return null;
  return ui.div(parts);
}

