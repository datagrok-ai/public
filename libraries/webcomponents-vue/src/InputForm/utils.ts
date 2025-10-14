import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {Observable} from 'rxjs';
import type {ValidationIcon} from '@datagrok-libraries/webcomponents/src/ValidationIcon/ValidationIcon';

//
// TODO: probably a separate type lib (?)
//
export type RestrictionType = 'disabled' | 'restricted' | 'info' | 'none';

export type ConsistencyInfo = {
  restriction: RestrictionType,
  inconsistent: boolean,
  assignedValue: any,
}

export interface ActionItem {
  actionName: string;
  action: string;
  additionalParams?: Record<string, any>;
}

export interface Advice {
  description: string;
  actions?: ActionItem[];
}

export type ValidationItem = string | Advice;

export interface ValidationResult {
  errors?: ValidationItem[];
  warnings?: ValidationItem[];
  notifications?: ValidationItem[];
}

export type ValidationIconInput = {
  validation?: ValidationResult,
  consistency?: ConsistencyInfo,
};

export interface SubscriptionLike {
  unsubscribe(): void;
}

export interface FuncCallInput<T = any> {
  root: HTMLElement;
  value: T | null | undefined;
  notify: boolean;
  enabled: boolean;
  onInput: ((cb: Function) => SubscriptionLike) | Observable<T>;
}


export const injectInputBaseStatus = (emit: Function, ioName: string, t: DG.InputBase) => {
  const icon = new (customElements.get('dg-validation-icon')!)() as ValidationIcon;
  icon.isScalar = DG.TYPES_SCALAR.has(t.property.type);
  icon.isDataFrame = t.property.type === DG.TYPE.DATA_FRAME;
  icon.addEventListener('consistency-reset', () => emit('consistencyReset', ioName));
  icon.addEventListener('action-request', (ev: any) => emit('actionRequested', ev.detail));
  icon.addEventListener('show-dataframe-diff', (ev: any) => showDFDiff(t.value, ev.detail));

  const wrapper = ui.element('i') as HTMLElement;
  wrapper.classList.add('rfv2-validation-icon');
  wrapper.appendChild(icon);
  t.addOptions(wrapper);

  function setStatus(status: {
    validation?: ValidationResult,
    consistency?: ConsistencyInfo,
  }) {
    const {validation} = status;

    icon.validationStatus = status;

    $(t.input).removeClass('d4-invalid d4-partially-invalid');
    if (validation?.errors && validation.errors.length)
      $(t.input).addClass('d4-invalid');
    else if (validation?.warnings && validation.warnings.length)
      $(t.input).addClass('d4-partially-invalid');
  }

  (t as any).setStatus = setStatus;
};

export interface FuncCallInputStatusable<T = any> extends FuncCallInput<T> {
  setStatus: (status: {
    validation?: ValidationResult,
    consistency?: ConsistencyInfo,
  }) => void;
}

export function isFuncCallInput<T = any>(arg: any): arg is FuncCallInput<T> {
  return arg && arg.root && arg.onInput;
}

export function isInputInjected(arg: any): arg is FuncCallInputStatusable {
  return arg?.setStatus && isFuncCallInput(arg);
}

function showDFDiff(df1?: DG.DataFrame, df2?: DG.DataFrame) {
  const idxName = '__compare_idx_col__'
  if (!df1 || !df2) {
    grok.shell.warning('One of dataframes is empty');
    return;
  }
  const cols1 = df1.columns.toList().map(col => col.name);
  const cols2 = df2.columns.toList().map(col => col.name);
  const df1c = df1.clone();
  df1c.columns.addNew(idxName, 'int').init((idx) => idx);
  const df2c = df2.clone();
  df2c.columns.addNew(idxName, 'int').init((idx) => idx);
  grok.data.compareTables(df1c, df2c, [idxName], [idxName], cols1, cols2, true);
}
