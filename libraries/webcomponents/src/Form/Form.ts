/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class Form extends HTMLElement {
  private formInst?: DG.InputForm;
  private currentSource?: DG.FuncCall;
  private initCalled = false;

  constructor() {
    super();
  }

  connectedCallback() {
    this.init();
  }

  disconnectedCallback() {
  }

  get funcCall() {
    return this.currentSource;
  }

  set funcCall(fc: DG.FuncCall | undefined) {
    this.currentSource = fc;
    if (this.formInst)
      this.formInst.source = this.currentSource;
    else
      this.init();
  }

  private async init() {
    if (this.initCalled)
      return;
    this.initCalled = true;
    // checking if a value has been updated during Form initialization
    const oldSource = this.currentSource;

    this.formInst = await DG.InputForm.forFuncCall(this.currentSource!, {twoWayBinding: true});

    if (oldSource !== this.currentSource)
      this.formInst.source = this.currentSource;

    this.attach();
  }

  private attach() {
    this.innerHTML = '';
    this.appendChild(this.formInst!.root);
  }
}

export interface FormT extends Form {};
