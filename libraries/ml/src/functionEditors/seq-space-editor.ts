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
      
      //@ts-ignore  
      this.funcParamsDiv = ui.inputs([
        this.tableInput,
        this.molColInput,
        this.methodInput,
        this.methodSettingsDiv,
        this.similarityMetricInput,
        this.plotEmbeddingsInput
      ], {style: {minWidth: '320px'}});
    }

  }
