import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IDimReductionParam, ITSNEOptions, IUMAPOptions, TSNEOptions, T_SNE, UMAP, UMAPOptions } from '../reduce-dimensionality';

export const SEQ_COL_NAMES = {
    [DG.SEMTYPE.MOLECULE]: 'Molecule',
    [DG.SEMTYPE.MACROMOLECULE]: 'Sequence'
}

export class SequenceSpaceBaseFuncEditor {
    tableInput: DG.InputBase;
    molColInput: DG.InputBase;
    methodInput: DG.InputBase;
    methodSettingsIcon: HTMLElement;
    methodSettingsDiv = ui.form([]);
    moleculesColDiv: HTMLDivElement = ui.div();
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
      this.molColInput.input.style.width = '160px';
      this.moleculesColDiv.append(this.molColInput.root);
  
      this.methodInput = ui.choiceInput('Method Name', UMAP, [UMAP, T_SNE], () => {
        this.createAlgorithmSettingsDiv(this.methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
      });
      this.methodInput.input.style.width = '160px';
      this.methodInput.captionLabel.style.width = '122px';
  
      this.methodSettingsIcon = ui.icons.settings(()=> {
        settingsOpened = !settingsOpened;
        if (!settingsOpened)
          ui.empty(this.methodSettingsDiv);
        else 
          this.createAlgorithmSettingsDiv(this.methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
      }, 'Modify methods parameters');
      this.methodSettingsIcon.style.lineHeight = '1.7';
      this.methodSettingsIcon.style.fontSize = '18px';
      this.methodSettingsDiv = ui.form([]);
      let settingsOpened = false;
    }
  
    createAlgorithmSettingsDiv(paramsForm: HTMLDivElement, params: UMAPOptions | TSNEOptions) {
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
        ui.empty(this.moleculesColDiv);
        this.molColInput = ui.columnInput(SEQ_COL_NAMES[semtype], this.tableInput.value!, this.tableInput.value!.columns.bySemType(semtype));
        this.molColInput.input.style.width = '160px';
        this.moleculesColDiv.append(this.molColInput.root);
    }
  }