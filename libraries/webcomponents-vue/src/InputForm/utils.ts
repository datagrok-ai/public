import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import type {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import $ from 'cash-dom';
import type {FuncCallInput} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';
import type {ValidationIcon} from '@datagrok-libraries/webcomponents/src/ValidationIcon/ValidationIcon';


export const injectInputBaseStatus = (emit: Function, ioName: string, t: DG.InputBase) => {
  const icon = new (customElements.get('dg-validation-icon')!)() as ValidationIcon;
  icon.isScalar = DG.TYPES_SCALAR.has(t.property.type);
  icon.addEventListener('consistency-reset', () => emit('consistencyReset', ioName));
  icon.addEventListener('action-request', (ev: any) => emit('actionRequested', ev.detail));

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
