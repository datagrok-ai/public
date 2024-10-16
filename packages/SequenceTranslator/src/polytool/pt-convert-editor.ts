import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {defaultErrorHandler} from '../utils/err-info';
import {PT_UI_DIALOG_CONVERSION, PT_UI_RULES_USED} from './const';
import {RuleInputs, RULES_PATH, RULES_STORAGE_NAME} from './pt-rules/pt-rules';

/** Inputs of polyToolConvert2 package function */
export enum P {
  table = 'table',
  seqCol = 'seqCol',
  generateHelm = 'generateHelm',
  chiralityEngine = 'chiralityEngine',
  rules = 'rules'
}

type PolyToolConvertInputs = {
  table: DG.InputBase<DG.DataFrame | null>;
  seqCol: DG.InputBase<DG.Column | null>;
  generateHelm: DG.InputBase<boolean>;
  chiralityEngine: DG.InputBase<boolean>;
  rules: { header: HTMLElement, form: HTMLDivElement };
}

export class PolyToolConvertFuncEditor {
  private inputs: PolyToolConvertInputs;
  private readonly ruleInputs: RuleInputs;

  protected constructor(
    private readonly call: DG.FuncCall,
  ) {
    this.ruleInputs = new RuleInputs(RULES_PATH, RULES_STORAGE_NAME, '.json');
  }

  private async initInputs(): Promise<void> {
    const getParam = (pName: string) => this.call.inputParams[pName];

    this.inputs = {
      table: (() => {
        const p = getParam(P.table);
        return ui.input.table(p.property.caption, {value: p.value});
      })(),
      seqCol: (() => {
        const p = getParam(P.seqCol);
        return ui.input.column(p.property.caption, {value: p.value, table: p.value.dataFrame});
      })(),
      generateHelm: ui.input.forProperty(getParam(P.generateHelm).property),
      chiralityEngine: ui.input.forProperty(getParam(P.chiralityEngine).property),
      rules: {
        header: ui.inlineText([PT_UI_RULES_USED]),
        form: await this.ruleInputs.getForm(),
      }
    };
  }

  static async create(call: DG.FuncCall): Promise<PolyToolConvertFuncEditor> {
    const editor = new PolyToolConvertFuncEditor(call);
    await editor.initInputs();
    return editor;
  }

  // -- Params --

  private async getParams(): Promise<{ [p: string]: any }> {
    return {
      table: this.inputs.table.value!,
      seqCol: this.inputs.seqCol.value!,
      generateHelm: this.inputs.generateHelm.value,
      chiralityEngine: this.inputs.chiralityEngine.value,
      rules: await this.ruleInputs.getActive(),
    };
  }

  public async showDialog(): Promise<DG.Column<string>> {
    const formDiv = ui.div([
      this.inputs.table,
      this.inputs.seqCol,
      this.inputs.generateHelm,
      this.inputs.chiralityEngine,
      this.inputs.rules.header,
      this.inputs.rules.form,
    ],
    {style: {minWidth: '320px'}});

    return new Promise((resolve, reject) => {
      ui.dialog({title: PT_UI_DIALOG_CONVERSION})
        .add(formDiv)
        .onOK(async () => {
          const callParams = await this.getParams();
          const res = (await this.call.func.prepare(callParams).call(true)).getOutputParamValue() as DG.Column<string>;
          resolve(res);
        })
        .onCancel(() => { reject(new Error('Cancelled by user')); })
        .show();
    });
  }

  // -- UI --

  public widget(): DG.Widget {
    throw new Error('not implemented');

    // const inputsForm = ui.inputs(Object.entries(this.inputs)
    //   .filter(([inputName, _input]) => !['table', 'sequence'].includes(inputName))
    //   .map(([_inputName, input]) => input));
    // const doBtn = ui.button('Convert', () => {
    //   (async () => {
    //     const callParams = this.getParams();
    //     await this.call.func.prepare(callParams).call(true);
    //   })()
    //     .catch((err: any) => { defaultErrorHandler(err, true); });
    // });
    // return DG.Widget.fromRoot(ui.divV([inputsForm, ui.div(doBtn)]));
  }
}
