import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {DimReductionEditorOptions} from '../functionEditors/dimensionality-reduction-editor';
import {MultiColumnDimReductionEditor}
  from '../multi-column-dimensionality-reduction/multi-column-dim-reduction-editor';


export class MCLEditor extends MultiColumnDimReductionEditor {
    public similarityThresholdInput: DG.InputBase<number | null>;
    public maxIterationsInput: DG.InputBase<number | null>;
    constructor(editorSettings: DimReductionEditorOptions = {}) {
      super(editorSettings);
      this.similarityThresholdInput = ui.intInput('Similarity threshold', 80);
      this.maxIterationsInput = ui.intInput('Max iterations', 5);
    }

    public getEditor(): HTMLElement {
      const div = ui.div([
        this.tableInput.root,
        this.columnsInputRoot,
        this.columnParamsEditorRoot,
        this.aggregationMethodInput.root,
        this.similarityThresholdInput.root,
        this.maxIterationsInput.root,
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
        weights: this.weightsEditor.weights,
        preprocessingFuncArgs: this.columnOptEditors.map((it) => it.preprocessingFunctionSettings),
        aggreaggregationMethod: this.aggregationMethodInput.value,
        threshold: this.similarityThresholdInput.value,
        maxIterations: this.maxIterationsInput.value ?? 5,
      };
    }
}
