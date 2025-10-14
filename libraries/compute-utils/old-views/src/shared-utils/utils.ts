import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {ACTIONS_COLUMN_NAME, AUTHOR_COLUMN_NAME, COMPLETE_COLUMN_NAME, DESC_COLUMN_NAME, EXP_COLUMN_NAME, FAVORITE_COLUMN_NAME, ID_COLUMN_NAME, STARTED_COLUMN_NAME, TAGS_COLUMN_NAME, TITLE_COLUMN_NAME} from '../../../shared-utils/consts';
import {FuncCallInput, isInputLockable, SubscriptionLike} from './input-wrappers';
import {ValidationResultBase, Validator, getValidationIcon, mergeValidationResults, nonNullValidator} from './validation';
import {FunctionView, RichFunctionView} from '../../';
import {getColumnName, getVisibleProps, setGridColumnsRendering} from '../../../shared-utils/history';
import {Observable} from 'rxjs';
import {getPropViewers} from '../../../shared-utils/utils';
import {SYNC_FIELD, SyncFields, syncParams, ValidationRequestPayload} from './consts';

export function isInputBase(input: FuncCallInput): input is DG.InputBase {
  const inputAny = input as any;
  return (inputAny.dart && DG.toJs(inputAny.dart) instanceof DG.InputBase);
}

export function getObservable<T>(onInput: (f: Function) => SubscriptionLike): Observable<T> {
  return new Observable((observer: any) => {
    const sub = onInput((val: T) => {
      observer.next(val);
    });
    return () => sub.unsubscribe();
  });
}

export const styleHistoryGrid = (
  grid: DG.Grid,
  isCompactMode: boolean,
  showInputsOnCards: boolean,
  showMetadataOnCards: boolean,
  func?: DG.Func,
  useOptions?: boolean,
) => {
  if (useOptions) {
    grid.setOptions({
      'showCurrentRowIndicator': true,
      'showCurrentCellOutline': false,
      'allowEdit': false,
      'allowBlockSelection': false,
      'showRowHeader': false,
      'showColumnLabels': !isCompactMode,
      'extendLastColumn': isCompactMode,
    });
  }

  grid.sort([STARTED_COLUMN_NAME], [false]);

  for (let i = 0; i < grid.columns.length; i++) {
    const col = grid.columns.byIndex(i);
    if (col && col.column?.type === DG.TYPE.DATE_TIME)
      col.format = 'MMM d, h:mm tt';
  }

  setGridColumnsRendering(grid);

  if (isCompactMode) {
    grid.columns.setVisible([ID_COLUMN_NAME]);

    grid.props.rowHeight = 70;
    grid.invalidate();
  } else {
    grid.props.rowHeight = 28;

    const tagCol = grid.dataFrame.getCol(TAGS_COLUMN_NAME);
    grid.columns.setVisible([
      EXP_COLUMN_NAME,
      FAVORITE_COLUMN_NAME,
      ACTIONS_COLUMN_NAME,
      ...showMetadataOnCards ? [STARTED_COLUMN_NAME]: [],
      ...showMetadataOnCards ? [AUTHOR_COLUMN_NAME]: [],
      ...showMetadataOnCards && tagCol.stats.missingValueCount < tagCol.length ? [TAGS_COLUMN_NAME]: [],
      ...showMetadataOnCards ? [TITLE_COLUMN_NAME]: [],
      ...showMetadataOnCards ? [DESC_COLUMN_NAME]: [],
      ...showInputsOnCards && func ? getVisibleProps(func)
        .map((key) => {
          const param = func.inputs.find((prop) => prop.name === key) ??
          func.outputs.find((prop) => prop.name === key);

          if (param)
            return param.caption ?? getColumnName(param.name);
          else
            return getColumnName(key);
        }): [],
    ]);
  }
};

export const styleHistoryFilters = (
  filters: DG.Viewer<DG.IFiltersSettings>,
  showMetadataColumns: boolean,
  showInputColumns: boolean,
  isHistory: boolean,
  func?: DG.Func,
) => {
  const currentDf = filters.dataFrame;
  const tagCol = currentDf.getCol(TAGS_COLUMN_NAME);

  const columnNames = [
    ...showMetadataColumns &&
    currentDf.getCol(EXP_COLUMN_NAME).categories.length > 1 ? [EXP_COLUMN_NAME]: [],
    ...showMetadataColumns &&
    (currentDf.col(FAVORITE_COLUMN_NAME)?.categories.length ?? 0) > 1 ?
      [FAVORITE_COLUMN_NAME]: [],
    ...showMetadataColumns ? [STARTED_COLUMN_NAME, COMPLETE_COLUMN_NAME]:[],
    ...showMetadataColumns &&
    currentDf.getCol(AUTHOR_COLUMN_NAME).categories.length > 1 ? [AUTHOR_COLUMN_NAME]: [],
    ...showMetadataColumns &&
    tagCol.stats.missingValueCount < tagCol.length ? [TAGS_COLUMN_NAME]: [],
    ...showMetadataColumns &&
    currentDf.getCol(TITLE_COLUMN_NAME).categories.length > 1 ? [TITLE_COLUMN_NAME]: [],
    ...showMetadataColumns &&
    currentDf.getCol(DESC_COLUMN_NAME).categories.length > 1 ? [DESC_COLUMN_NAME]: [],
    ...func && showInputColumns ? getVisibleProps(func)
      .map((key) => {
        const param = func.inputs.find((prop) => prop.name === key) ??
      func.outputs.find((prop) => prop.name === key);

        if (param)
          return param.caption ?? getColumnName(param.name);
        else
          return getColumnName(key);
      })
      .filter((columnName) => {
        return !isHistory ||
        currentDf.getCol(columnName).categories.length > 1;
      })
      .map((columnName) => columnName): [],
  ];
  if (columnNames.length > 0) {
    ui.setDisplay(filters.root, true);
    filters.setOptions({columnNames, 'showHeader': false, 'showBoolCombinedFilter': true});
  } else
    ui.setDisplay(filters.root, false);


  return columnNames.length > 0;
};

export const showHelpWithDelay = async (helpContent: string) => {
  grok.shell.windows.help.visible = true;
  // Workaround to deal with help panel bug
  await new Promise((resolve) => setTimeout(resolve, 100));
  grok.shell.windows.help.showHelp(ui.markdown(helpContent));
};

export const categoryToDfParamMap = (func: DG.Func) => {
  const map = {
    inputs: {} as Record<string, DG.Property[]>,
    outputs: {} as Record<string, DG.Property[]>,
  };

  func.inputs
    .filter((inputProp) =>
      inputProp.propertyType === DG.TYPE.DATA_FRAME &&
      getPropViewers(inputProp).config.length !== 0,
    )
    .forEach((p) => {
      const category = p.category === 'Misc' ? 'Input': p.category;

      if (map.inputs[category])
        map.inputs[category].push(p);
      else
        map.inputs[category] = [p];
    });

  func.outputs
    .forEach((p) => {
      const category = p.category === 'Misc' ? 'Output': p.category;

      if (p.propertyType === DG.TYPE.DATA_FRAME &&
        getPropViewers(p).config.length === 0) return;

      if (map.outputs[category])
        map.outputs[category].push(p);
      else
        map.outputs[category] = [p];
    });

  return map;
};

export const createPartialCopy = async (call: DG.FuncCall) => {
  const previousId = call.id;
  // grok.functions.eval creates an ID.
  // So we should control it and overwrite ID by null again if necessary.

  const callCopy: DG.FuncCall = DG.Func.byName(call.func.nqName)
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

export const getValidators = async (funcCall: DG.FuncCall, isInput: SyncFields = SYNC_FIELD.INPUTS) => {
  const params = [...funcCall[syncParams[isInput]].values()];
  const resolvedValidators = await Promise.all(
    params
      .filter((param) => !!param.property.options.validatorFunc)
      .map(async (param) => {
        const func = DG.Func.byName(param.property.options.validatorFunc);
        const call = func.prepare({params: JSON.parse(param.property.options.validatorFuncOptions || '{}')});
        await call.call();

        return [param.name, call.outputs.validator as Validator] as const;
      }));

  return resolvedValidators.reduce((acc, [name, validator]) => {
    acc[name] = validator;
    return acc;
  }, {} as Record<string, Validator>);
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

export const validate = async (
  payload: ValidationRequestPayload,
  paramNames: string[],
  signal: AbortSignal,
  isInput: SyncFields,
  context: {view?: RichFunctionView, funcCall: DG.FuncCall, lastCall?: DG.FuncCall},
  validarors: Record<string, Validator>,
) => {
  const {view, funcCall, lastCall} = context;

  const validationItems = await Promise.all(paramNames.map(async (name) => {
    const v = isInput === SYNC_FIELD.INPUTS ? funcCall.inputs[name]: funcCall.outputs[name];
    // not allowing null anywhere
    const standardMsgs = await nonNullValidator(v, {
      param: name,
      funcCall: funcCall,
      lastCall: lastCall,
      signal,
      isNewOutput: !!payload.isNewOutput,
      isRevalidation: payload.isRevalidation,
      view: view!,
    });
    let customMsgs;
    const customValidator = validarors[name];
    if (customValidator) {
      customMsgs = await customValidator(v, {
        param: name,
        funcCall: funcCall,
        lastCall: lastCall,
        signal,
        isNewOutput: !!payload.isNewOutput,
        isRevalidation: payload.isRevalidation,
        context: payload.context,
        view: view!,
      });
    }
    // output params could not be nulls, DG will complain
    const isNullable = isInput === SYNC_FIELD.INPUTS && funcCall.inputParams[name].property.options.nullable;
    return [name, mergeValidationResults(
      ...isNullable ? []: [standardMsgs],
      customMsgs,
    )] as const;
  }));
  return Object.fromEntries(validationItems);
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

export const injectLockIcons = (
  t: FuncCallInput,
  onUnlock: Function,
  onUndo: Function,
  onWarningClick: Function,
) => {
  t.root.addEventListener('click', () => onUnlock());

  const lockIcon = ui.iconFA('lock');
  $(lockIcon).addClass('rfv-icon-lock');
  $(lockIcon).css({color: `var(--grey-2)`});

  const unlockIcon = ui.iconFA('lock-open');
  $(unlockIcon).addClass('rfv-icon-unlock');
  $(unlockIcon).css({color: `var(--grey-2)`});

  const resetIcon = ui.iconFA('undo', (e: MouseEvent) => onUndo(e), 'Reset value to computed value');
  $(resetIcon).addClass('rfv-icon-undo');
  $(resetIcon).css({color: `var(--blue-2)`});

  const warningIcon = ui.iconFA('exclamation-circle', null);
  ui.tooltip.bind(warningIcon, () => onWarningClick());
  $(warningIcon).addClass('rfv-icon-warning');
  $(warningIcon).css({color: `var(--orange-2)`});

  function defaultPlaceLockStateIcons(
    lockIcon: HTMLElement,
    unlockIcon: HTMLElement,
    resetIcon: HTMLElement,
    warningIcon: HTMLElement,
  ) {
    // If custom input is not DG.InputBase instance then do nothing
    if (!isInputBase(t)) return;

    t.addOptions(lockIcon);
    t.addOptions(unlockIcon);
    t.addOptions(resetIcon);
    t.addOptions(warningIcon);
  }

  const tAny = (t as any);
  // if no custom place for lock state icons is provided then use default placing
  if (!tAny.placeLockStateIcons)
    tAny.placeLockStateIcons = defaultPlaceLockStateIcons;
  tAny.placeLockStateIcons(lockIcon, unlockIcon, resetIcon, warningIcon);
};
