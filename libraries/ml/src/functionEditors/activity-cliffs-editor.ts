import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IDimReductionParam, TSNEOptions, UMAPOptions } from '../reduce-dimensionality';
import { SequenceSpaceBaseFuncEditor } from './seq-space-base-editor';

export class ActivityCliffsFunctionEditor extends SequenceSpaceBaseFuncEditor {

    activitiesInput: DG.InputBase;
    similarityInput: DG.InputBase;
    funcParamsDiv: HTMLDivElement;
    activitiesColDiv: HTMLDivElement = ui.div();

    get funcParams(): any {
      return {table: this.tableInput.value!, molecules: this.molColInput.value!, activities: this.activitiesInput.value!,
        similarity: this.similarityInput.value!, methodName: this.methodInput.value!, options: this.algorithmOptions};
    }
  
    get paramsUI(): HTMLDivElement{
      return this.funcParamsDiv;
    }
  
    constructor(semtype: DG.SemType){
      super(semtype);
      this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, this.tableInput.value!.columns.byIndex(0));
      this.activitiesInput.input.style.width = '160px';
      this.activitiesColDiv.append(this.activitiesInput.root);
      
      this.similarityInput = ui.intInput('Similarity', 80);
      
      this.funcParamsDiv = ui.form([
        this.tableInput,
        //@ts-ignore
        this.moleculesColDiv,
        //@ts-ignore
        this.activitiesColDiv,
        this.similarityInput,
        //@ts-ignore
        ui.divH([this.methodSettingsIcon, this.methodInput]),
        //@ts-ignore
        this.methodSettingsDiv
      ])
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
        super.onTableInputChanged(semtype);
        ui.empty(this.activitiesColDiv);
        this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, this.tableInput.value!.columns.byIndex(0));
        this.activitiesInput.input.style.width = '160px';
        this.activitiesColDiv.append(this.activitiesInput.root);
    }
  }