import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import wu from 'wu';
import type ExcelJS from 'exceljs';
import type html2canvas from 'html2canvas';
import {AUTHOR_COLUMN_NAME, VIEWER_PATH, viewerTypesMapping} from './consts';
import {FuncCallInput, isInputLockable} from './input-wrappers';
import {ValidationResultBase, getValidationIcon} from './validation';
import {FunctionView, RichFunctionView} from '../function-views';

export const updateIndicatorWithText = (element: HTMLElement, updating: boolean, text?: string) => {
  ui.setUpdateIndicator(element, updating);
  const updatingLabel = element.querySelector('.d4-update-shadow .ui-label');
  if (updating && text && updatingLabel) {
    updatingLabel.textContent = text;
  }
};

export const createPartialCopy = async (call: DG.FuncCall) => {
  const previousId = call.id;
  // grok.functions.eval creates an ID.
  // So we should control it and overwrite ID by null again if necessary.

  const callCopy: DG.FuncCall = (await grok.functions.eval(call.func.nqName))
    //@ts-ignore
    .prepare([...call.inputs].reduce((acc, [key, val]) => {
      acc[key] = val;
      return acc;
    }, {} as Record<string, any>));
  call.options.forEach((key: string) => callCopy.options[key] = call.options[key]);

  //@ts-ignore
  if (!previousId) callCopy.id = null;

  return callCopy;
};

export const isIncomplete = (run: DG.FuncCall) => {
  return !getStartedOrNull(run) || !run.id;
};

export const getStartedOrNull = (run: DG.FuncCall) => {
  try {
    return run.started;
  } catch {
    return null;
  }
};

export const extractStringValue = (run: DG.FuncCall, key: string) => {
  if (key === AUTHOR_COLUMN_NAME) return run.author?.friendlyName ?? grok.shell.user.friendlyName;

  const val =
  (run as any)[key] ??
  run.inputs[key] ??
  run.outputs[key] ??
  run.options[key] ??
  null;

  return val?.toString() ?? '';
};

export const getMainParams = (func: DG.Func): string[] | null => {
  return func.options['mainParams'] ? JSON.parse(func.options['mainParams']): null;
};

export const camel2title = (camelCase: string) => camelCase
  .replace(/([A-Z])/g, (match) => ` ${match.toLowerCase()}`)
  .trim()
  .replace(/^./, (match) => match.toUpperCase());

export function isInputBase(input: FuncCallInput): input is DG.InputBase {
  const inputAny = input as any;
  return (inputAny.dart && DG.toJs(inputAny.dart) instanceof DG.InputBase);
}

export const deepCopy = (call: DG.FuncCall) => {
  const previousId = call.id;
  // FuncCall.clone() creates an ID for original (!) call if it was null.
  // So we should control it and overwrite ID by null again if necessary.
  const deepClone = call.clone();

  //@ts-ignore
  if (!previousId) deepClone.id = null;

  call.options.forEach((key: string) => deepClone.options[key] = call.options[key]);

  const dfOutputs = wu(call.outputParams.values())
    .filter((output) =>
      output.property.propertyType === DG.TYPE.DATA_FRAME &&
      !!call.outputs[output.name],
    );
  for (const output of dfOutputs)
    deepClone.outputs[output.name] = call.outputs[output.name].clone();

  const dfInputs = wu(call.inputParams.values())
    .filter((input) =>
      input.property.propertyType === DG.TYPE.DATA_FRAME &&
      !!call.inputs[input.name],
    );
  for (const input of dfInputs)
    deepClone.inputs[input.name] = call.inputs[input.name].clone();

  return deepClone;
};

export const getPropViewers = (prop: DG.Property): {name: string, config: Record<string, string | boolean>[]} => {
  const viewersRawConfig = prop.options[VIEWER_PATH];
  return viewersRawConfig ?
  // true and false values are retrieved as string, so we parse them separately
    {name: prop.name, config: JSON.parse(viewersRawConfig, (k, v) => {
      if (v === 'true') return true;
      if (v === 'false') return false;
      // Converting internal Dart labels to JS DG.VIEWER labels
      if (k === 'type') return viewerTypesMapping[v] || v;

      if (!k.toLowerCase().includes('color')) {
        const parsed = Number.parseFloat(v);

        if (!Number.isNaN(parsed))
          return parsed;
      }

      return v;
    })}:
    {name: prop.name, config: []};
};

export const getFuncRunLabel = (func: DG.Func) => {
  return func.options['runLabel'];
};

export const injectLockStates = (input: FuncCallInput) => {
  // if custom lock state methods are available then use them
  if (isInputLockable(input)) return;

  function setDisabledDefault() {
    input.enabled = false;
    $(input.root).removeClass('rfv-restricted-unlocked-input');
    $(input.root).removeClass('rfv-inconsistent-input');
    $(input.root).removeClass('rfv-restricted-input');
  }

  function setRestrictedDefault() {
    input.enabled = false;
    if (isInputBase(input)) (input.input as HTMLInputElement).disabled = false; ;
    $(input.root).addClass('rfv-restricted-input');
    $(input.root).removeClass('rfv-restricted-unlocked-input');
    $(input.root).removeClass('rfv-inconsistent-input');
  }

  function setRestrictedUnlockedDefault() {
    input.enabled = true;
    $(input.root).addClass('rfv-restricted-unlocked-input');
    $(input.root).removeClass('rfv-restricted-input');
    $(input.root).removeClass('rfv-inconsistent-input');
  }

  function setInconsistentDefault() {
    input.enabled = true;
    $(input.root).addClass('rfv-inconsistent-input');
    $(input.root).removeClass('rfv-restricted-input');
    $(input.root).removeClass('rfv-restricted-unlocked-input');
  }

  function setUserInputDefault() {
    input.enabled = true;
    $(input.root).removeClass('rfv-inconsistent-input');
    $(input.root).removeClass('rfv-restricted-input');
    $(input.root).removeClass('rfv-restricted-unlocked-input');
  }

  const inputAny = input as any;
  inputAny.setDisabled = setDisabledDefault;
  inputAny.setRestricted = setRestrictedDefault;
  inputAny.setRestrictedUnlocked = setRestrictedUnlockedDefault;
  inputAny.setInconsistent = setInconsistentDefault;
  inputAny.setUserInput = setUserInputDefault;
};

export const inputBaseAdditionalRenderHandler = (val: DG.FuncCallParam, t: DG.InputBase) => {
  const prop = val.property;

  $(t.root).css({
    'width': `calc(${prop.options['block'] ?? '100'}% - ${prop.options['block'] ? '2': '0'}px)`,
    'box-sizing': 'border-box',
  });
};

export const updateOutputValidationSign = (
  sign: readonly [HTMLElement, HTMLElement],
  messages: ValidationResultBase | undefined,
):readonly [HTMLElement, HTMLElement] => {
  const newSign = getValidationIcon(messages);
  sign[0].replaceWith(newSign[0]);
  sign[1].replaceWith(newSign[1]);

  return newSign;
};

export const injectInputBaseValidation = (t: DG.InputBase) => {
  const validationIndicator = ui.element('i');
  t.addOptions(validationIndicator);
  function setValidation(messages: ValidationResultBase | undefined) {
    while (validationIndicator.firstChild && validationIndicator.removeChild(validationIndicator.firstChild));
    const [icon, popover] = getValidationIcon(messages);
    if (icon && popover) {
      validationIndicator.appendChild(icon);
      validationIndicator.appendChild(popover);
    }

    t.input.classList.remove('d4-invalid');
    t.input.classList.remove('d4-partially-invalid');
    if (
      (messages?.errors && messages.errors.length) ||
      (messages?.warnings && messages.warnings.length) ||
      (messages?.notifications && messages.notifications.length) ||
      messages?.pending
    )
      $(validationIndicator).css('display', 'flex');
    else
      $(validationIndicator).hide();

    if (messages?.errors && messages.errors.length)
      t.input.classList.add('d4-invalid');
    else if (messages?.warnings && messages.warnings.length)
      t.input.classList.add('d4-partially-invalid');
  }
  (t as any).setValidation = setValidation;
};

export const scalarsToSheet =
  (sheet: ExcelJS.Worksheet, scalars: { caption: string, value: string, units: string }[]) => {
    sheet.addRow(['Parameter', 'Value', 'Units']).font = {bold: true};
    scalars.forEach((scalar) => {
      sheet.addRow([scalar.caption, scalar.value, scalar.units]);
    });

    sheet.getColumn(1).width = Math.max(
      ...scalars.map((scalar) => scalar.caption.toString().length), 'Parameter'.length,
    ) * 1.2;
    sheet.getColumn(2).width = Math.max(
      ...scalars.map((scalar) => scalar.value.toString().length), 'Value'.length) * 1.2;
    sheet.getColumn(3).width = Math.max(
      ...scalars.map((scalar) => scalar.units.toString().length), 'Units'.length) * 1.2;
  };

let dfCounter = 0;
export const dfToSheet = (sheet: ExcelJS.Worksheet, df: DG.DataFrame, column?: number, row?: number) => {
  const columnKey = sheet.getColumn(column ?? 1).letter;
  const tableConfig = {
    name: `ID_${dfCounter.toString()}`,
    ref: `${columnKey}${row ?? 1}`,
    columns: df.columns.toList().map((col) => ({name: col.name, filterButton: false})),
    rows: new Array(df.rowCount).fill(0).map((_, idx) => [...df.row(idx).cells].map((cell) => cell.value)),
  };
  sheet.addTable(tableConfig);
  sheet.columns.forEach((col) => {
    col.width = 25;
    col.alignment = {wrapText: true};
  });
  dfCounter++;
};

// additional JSON converions, view is need for files
export async function fcToSerializable(fc: DG.FuncCall, view: FunctionView | RichFunctionView) {
  const inputs: Record<string, any> = {};
  for (const [name, value] of Object.entries(fc.inputs)) {
    const {property} = view.funcCall.inputParams[name];
    inputs[name] = await fcInputToSerializable(property, value, view);
  }
  return {
    inputs,
    outputs: fc.outputs,
  };
}

async function fcInputToSerializable(property: DG.Property, value: any, view: FunctionView | RichFunctionView) {
  if ((property.propertyType as any) === 'file' && (view as any)!.getInput) {
    const fileInput = (view as any).getInput(property.name);
    return fileInput.value.arrayBuffer();
  }
  return value;
}

export async function fcInputFromSerializable(propertyType: string, value: any) {
  if (propertyType === 'file')
    return new File([value], '');

  return value;
}
