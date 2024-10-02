import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {DimReductionEditorOptions} from '../functionEditors/dimensionality-reduction-editor';
import {MultiColumnDimReductionEditor}
  from '../multi-column-dimensionality-reduction/multi-column-dim-reduction-editor';
import {getGPUAdapterDescription} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {defaultMCLOptions} from './markov-cluster';


export class MCLEditor extends MultiColumnDimReductionEditor {
    public similarityThresholdInput: DG.InputBase<number | null>;
    public maxIterationsInput: DG.InputBase<number | null>;
    public useWebGPUInput: DG.InputBase<boolean | null>;
    public inflateInput: DG.InputBase<number | null>;
    public minClusterSizeInput: DG.InputBase<number | null>;
    constructor(editorSettings: DimReductionEditorOptions = {}) {
      super(editorSettings);
      this.similarityThresholdInput = ui.input.int('Similarity Threshold', {value: 80});
      this.maxIterationsInput = ui.input.int('Max Iterations', {value: 5});
      this.useWebGPUInput = ui.input.bool('Use WebGPU', {value: false});
      this.inflateInput = ui.input.float('Inflation Factor', {value: defaultMCLOptions.inflateFactor});
      this.minClusterSizeInput = ui.input.int('Min Cluster Size', {value: 5});
      getGPUAdapterDescription().then((desc) => {
        if (desc) {
          this.useWebGPUInput.setTooltip(`Use webGPU for MCL calculation (${desc})`);
          this.useWebGPUInput.value = true;
        } else {
          this.useWebGPUInput.value = false;
          this.useWebGPUInput.setTooltip('WebGPU is not available');
          this.useWebGPUInput.enabled = false;
        }
      });
    }

    public getEditor(): HTMLElement {
      const div = ui.div([
        this.tableInput.root,
        this.columnsInputRoot,
        this.columnParamsEditorRoot,
        this.aggregationMethodInput.root,
        this.similarityThresholdInput.root,
        this.inflateInput.root,
        this.maxIterationsInput.root,
        this.minClusterSizeInput.root,
        this.useWebGPUInput.root,
      ], {style: {minWidth: '420px'}, classes: 'ui-form'});
      return div;
    }

    public get params() {
      return {
        table: this.tableInput.value!,
        columns: this.columnsInput.value!,
        methodName: this.methodInput.value!,
        preprocessingFunctions: this.columnOptEditors.map((it) => it.preProcessingFunction),
        distanceMetrics: this.columnOptEditors.map((it) => it.similarityMetricInput.value!),
        weights: this.columnOptEditors.map((it) => it.weight!),
        preprocessingFuncArgs: this.columnOptEditors.map((it) => it.preprocessingFunctionSettings),
        aggreaggregationMethod: this.aggregationMethodInput.value,
        threshold: this.similarityThresholdInput.value,
        maxIterations: this.maxIterationsInput.value ?? 5,
        useWebGPU: this.useWebGPUInput.value ?? false,
        inflateFactor: this.inflateInput.value ?? defaultMCLOptions.inflateFactor,
        minClusterSize: this.minClusterSizeInput.value ?? 5,
      };
    }
}
