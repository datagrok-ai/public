import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {VIEWER_PATH, viewerTypesMapping} from './consts';
import $ from 'cash-dom';
import wu from 'wu';
import {FuncCallInput, isInputLockable} from './input-wrappers';
import {ValidationResultBase, getValidationIcon} from './validation';

export const deepCopy = (call: DG.FuncCall) => {
  const deepClone = call.clone();

  const dfOutputs = wu(call.outputParams.values() as DG.FuncCallParam[])
    .filter((output) => output.property.propertyType === DG.TYPE.DATA_FRAME);
  for (const output of dfOutputs)
    deepClone.outputs[output.name] = call.outputs[output.name].clone();

  const dfInputs = wu(call.inputParams.values() as DG.FuncCallParam[])
    .filter((input) => input.property.propertyType === DG.TYPE.DATA_FRAME);
  for (const input of dfInputs)
    deepClone.inputs[input.name] = call.inputs[input.name].clone();

  return deepClone;
};

export const boundImportFunction = (func: DG.Func): string | undefined => {
  return func.options['getRealData'];
};

export const getPropViewers = (prop: DG.Property): {name: string, config: Record<string, string | boolean>[]} => {
  const viewersRawConfig = prop.options[VIEWER_PATH];
  return (viewersRawConfig !== undefined) ?
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
    $(input.root).removeClass('rfv-inconsistent-input rfv-inconsistent-warning-input');
    $(input.root).removeClass('rfv-restricted-input');
  }

  function setRestrictedDefault() {
    input.enabled = false;
    $(input.root).removeClass('rfv-inconsistent-input rfv-inconsistent-warning-input');
    $(input.root).addClass('rfv-restricted-input');
  }

  function setInconsistentWarnDefault() {
    input.enabled = true;
    $(input.root).addClass('rfv-inconsistent-warning-input');

    $(input.root).addClass('rfv-inconsistent-input');
    $(input.root).removeClass('rfv-restricted-input');
  }

  function setInconsistentDefault() {
    input.enabled = true;
    $(input.root).removeClass('rfv-inconsistent-warning-input');

    $(input.root).addClass('rfv-inconsistent-input');
    $(input.root).removeClass('rfv-restricted-input');
  }

  const inputAny = input as any;
  inputAny.setDisabled = setDisabledDefault;
  inputAny.setRestricted = setRestrictedDefault;
  inputAny.setInconsistentWarn = setInconsistentWarnDefault;
  inputAny.setInconsistent = setInconsistentDefault;
};

export const inputBaseAdditionalRenderHandler = (val: DG.FuncCallParam, t: DG.InputBase) => {
  const prop = val.property;

  $(t.root).css({
    'width': `${prop.options['block'] ?? '100'}%`,
    'box-sizing': 'border-box',
  });
    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13004
    t.captionLabel.firstChild!.replaceWith(ui.span([prop.caption ?? prop.name]));
    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13005
    if (prop.options['units']) t.addPostfix(prop.options['units']);
};

export const injectInputBaseValidation = (t: DG.InputBase) => {
  const validationIndicator = ui.div('', {style: {display: 'flex'}});
  t.addOptions(validationIndicator);
  function setValidation(messages: ValidationResultBase | undefined) {
    while (validationIndicator.firstChild && validationIndicator.removeChild(validationIndicator.firstChild));
    const icon = getValidationIcon(messages);
    if (icon)
      validationIndicator.appendChild(icon);

    t.input.classList.remove('d4-invalid');
    t.input.classList.remove('d4-partially-invalid');
    if (messages?.errors)
      t.input.classList.add('d4-invalid');
    else if (messages?.warnings)
      t.input.classList.add('d4-partially-invalid');
  }
  (t as any).setValidation = setValidation;
};
