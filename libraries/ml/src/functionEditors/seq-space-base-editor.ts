import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IDimReductionParam, ITSNEOptions, IUMAPOptions, TSNEOptions, T_SNE, UMAP, UMAPOptions } from '../reduce-dimensionality';

export const SEQ_COL_NAMES = {
    [DG.SEMTYPE.MOLECULE]: 'Molecules',
    [DG.SEMTYPE.MACROMOLECULE]: 'Sequences'
}

export class SequenceSpaceBaseFuncEditor {
    tableInput: DG.InputBase;
    molColInput: DG.InputBase;
    molColInputRoot: HTMLElement;
    methodInput: DG.InputBase;
    methodSettingsIcon: HTMLElement;
    methodSettingsDiv = ui.inputs([]);
    methodsParams: {[key: string]: UMAPOptions | TSNEOptions} = {
      [UMAP]: new UMAPOptions(),
      [T_SNE]: new TSNEOptions()
    };
  
    get algorithmOptions(): IUMAPOptions | ITSNEOptions {
      const algorithmParams: UMAPOptions | TSNEOptions = this.methodsParams[this.methodInput.value!];
      const options: any = {};
      Object.keys(algorithmParams).forEach((key: string) => {
        if ((algorithmParams as any)[key].value != null) 
          options[key] = (algorithmParams as any)[key].value;
      });
      return options;
    }
  
    constructor(semtype: DG.SemType){
      this.tableInput = ui.tableInput('Table', grok.shell.tv.dataFrame, undefined, () => {
        this.onTableInputChanged(semtype);
      });
  
      this.molColInput = ui.columnInput(SEQ_COL_NAMES[semtype], this.tableInput.value!, this.tableInput.value!.columns.bySemType(semtype));
      this.molColInputRoot = this.molColInput.root;
      this.methodInput = ui.choiceInput('Method', UMAP, [UMAP, T_SNE], () => {
        if(settingsOpened) {
            this.createAlgorithmSettingsDiv(this.methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
        }
      });
  
      this.methodSettingsIcon = ui.icons.settings(()=> {
        settingsOpened = !settingsOpened;
        if (!settingsOpened)
          ui.empty(this.methodSettingsDiv);
        else 
          this.createAlgorithmSettingsDiv(this.methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
      }, 'Modify methods parameters');
      this.methodInput.root.classList.add('ml-dim-reduction-settings-input');
      this.methodInput.root.prepend(this.methodSettingsIcon);
      this.methodSettingsDiv = ui.inputs([]);
      let settingsOpened = false;
    }
  
    createAlgorithmSettingsDiv(paramsForm: HTMLDivElement, params: UMAPOptions | TSNEOptions): HTMLElement {
      ui.empty(paramsForm);
      Object.keys(params).forEach((it: any) => {
        const param: IDimReductionParam = (params as any)[it];
        const input = ui.floatInput(param.uiName, param.value, () => {
          param.value = input.value;
        });
        ui.tooltip.bind(input.root, param.tooltip);
        paramsForm.append(input.root);
      });
      return paramsForm;
    }

    onTableInputChanged(semtype: DG.SemType) {
        this.molColInput = ui.columnInput(SEQ_COL_NAMES[semtype], this.tableInput.value!, this.tableInput.value!.columns.bySemType(semtype));
        ui.empty(this.molColInputRoot);
        Array.from(this.molColInput.root.children).forEach((it) => this.molColInputRoot.append(it));
    }
  }
