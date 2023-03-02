import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IDimReductionParam, TSNEOptions, UMAPOptions } from '../reduce-dimensionality';
import { SequenceSpaceBaseFuncEditor } from './seq-space-base-editor';
import { SEQ_SPACE_SIMILARITY_METRICS } from '../distance-metrics-methods';

export class SequenceSpaceFunctionEditor extends SequenceSpaceBaseFuncEditor {

    similarityMetricInput: DG.InputBase;
    plotEmbeddingsInput: DG.InputBase;
    funcParamsDiv: HTMLDivElement;
  
    get funcParams(): any {
      return {table: this.tableInput.value!, molecules: this.molColInput.value!, methodName: this.methodInput.value!,
        similarityMetric: this.similarityMetricInput.value!, plotEmbeddings: this.plotEmbeddingsInput.value!, 
        options: this.algorithmOptions};
    }
  
    get paramsUI(): HTMLDivElement{
      return this.funcParamsDiv;
    }
  
    constructor(semtype: DG.SemType){
      super(semtype);
      this.similarityMetricInput = ui.choiceInput('Similarity metric', 'Tanimoto', SEQ_SPACE_SIMILARITY_METRICS);
      
      this.plotEmbeddingsInput = ui.boolInput('Plot Embeddings', true);
      this.plotEmbeddingsInput.captionLabel.style.width = '130px';
      this.plotEmbeddingsInput.input.style.marginLeft = '0px';
      
      this.funcParamsDiv = ui.form([
        this.tableInput,
        //@ts-ignore
        this.moleculesColDiv,
        //@ts-ignore
        ui.divH([this.methodSettingsIcon, this.methodInput]),
        //@ts-ignore
        this.methodSettingsDiv,
        this.similarityMetricInput,
        this.plotEmbeddingsInput
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
  }