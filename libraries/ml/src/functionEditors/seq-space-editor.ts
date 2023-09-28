import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ITSNEOptions, IUMAPOptions } from '../reduce-dimensionality';
import { SequenceSpaceBaseFuncEditor } from './seq-space-base-editor';

export interface ISimilaritySpaceParams {
  table: DG.DataFrame;
  molecules: DG.Column;
  methodName: string;
  similarityMetric: string;
  plotEmbeddings: boolean;
  sparseMatrixThreshold?: number;
  options: IUMAPOptions | ITSNEOptions;
}

export class SequenceSpaceFunctionEditor extends SequenceSpaceBaseFuncEditor {

  plotEmbeddingsInput: DG.InputBase;
  funcParamsDiv: HTMLDivElement;

  get funcParams(): ISimilaritySpaceParams {
    return {
      table: this.tableInput.value!,
      molecules: this.molColInput.value!,
      methodName: this.methodInput.value!,
      similarityMetric: this.similarityMetricInput.value!,
      sparseMatrixThreshold: this.similarityThresholdInput.value ?? 0.8,
      plotEmbeddings: this.plotEmbeddingsInput.value!,
      options: this.algorithmOptions
    };
  }

  get paramsUI(): HTMLDivElement{
    return this.funcParamsDiv;
  }

  constructor(semtype: DG.SemType){
    super(semtype);

    this.plotEmbeddingsInput = ui.boolInput('Plot Embeddings', true);

    //@ts-ignore
    this.funcParamsDiv = ui.inputs([
      this.tableInput,
      this.molColInput,
      this.methodInput,
      this.similarityThresholdInput,
      this.methodSettingsDiv,
      this.similarityMetricInput,
      this.plotEmbeddingsInput
    ], {style: {minWidth: '320px'}});
  }
}
