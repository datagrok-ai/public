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
    if (!fc) {
      ui.empty(this);
      return;
    }

    if (!this.funcCall) {
      this.replaceFunc(fc);
      return;
    }

    // if (this.funcCall.func.id === fc.func.id) {
    //   if (this.funcCall.id !== fc.id) this.formInst!.source = this.funcCall;
    // } else
      this.replaceFunc(fc);
  }

  private async replaceFunc(funcCall: DG.FuncCall) {
    const fcInputs = [...funcCall.inputs.entries()] as [string, any][];
    this.formInst = await DG.InputForm.forFuncCall(funcCall, {twoWayBinding: false});
    this.attach();

    setTimeout(() => {
      fcInputs.forEach(([paramName, paramValue]) => {
        const input = this.formInst!.getInput(paramName);
        input.notify = false;
        try {
          input.value = paramValue;
        } finally {
          input.notify = true;
        }
      });
      this.formInst!.onInputChanged
        .subscribe((event) => this.dispatchEvent(new CustomEvent('input-changed', {detail: event})));
    }, 0);
  }

  private attach() {
    ui.empty(this);
    if (this.formInst)
      this.appendChild(this.formInst.root);
    this.dispatchEvent(new CustomEvent('form-replaced', {detail: this.formInst}));
  }
}

export interface InputFormT extends InputForm {};
