/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class InputForm extends HTMLElement {
  private formInst?: DG.InputForm;

  constructor() {
    super();
  }

  connectedCallback() {
  }

  disconnectedCallback() {
  }

  get funcCall() {
    return this.formInst?.source;
  }

  set funcCall(fc: DG.FuncCall | undefined) {
    ui.empty(this);
    this.replaceFunc(fc);
  }

  private async replaceFunc(funcCall?: DG.FuncCall) {
    if (!funcCall)
      this.formInst = undefined;
    else {
      this.formInst = await DG.InputForm.forFuncCall(funcCall, {twoWayBinding: true, skipDefaultInit: true} as any);
      this.appendChild(this.formInst.root);
    }
    this.dispatchEvent(new CustomEvent('form-replaced', {detail: this.formInst}));
  }
}

export interface InputFormT extends InputForm {};
