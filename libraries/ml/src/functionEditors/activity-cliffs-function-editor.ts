import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ColumnInputOptions} from '@datagrok-libraries/utils/src/type-declarations';
import {DimReductionBaseEditor, DimReductionEditorOptions, DimReductionParams} from './dimensionality-reduction-editor';
import {IDimReductionParam} from '../multi-column-dimensionality-reduction/multi-column-dim-reducer';
import {MCLMethodName} from '../MCL';

export type ActivityCliffsParams = DimReductionParams & {
    activities: DG.Column;
    similarityThreshold: number;
}

export class MCLOptions {
    maxIterations: IDimReductionParam =
    {uiName: 'Max iterations', value: 0, tooltip: `Maximum iterations for MCL process.Default is 
      0 which will construct the clusters with plain sparse matrix. Values greater than 0 will 
      perform MCL with the given number of iterations and will result in trans-cluster activity cliff lines.`};
    useWebGPU: IDimReductionParam<boolean> = {
      uiName: 'Use WebGPU', value: false, tooltip: 'Use WebGPU for MCL calculation. (Experimental)',
      type: 'boolean',
    }
    constructor() {}
}

export class ActivityCliffsEditor extends DimReductionBaseEditor {
    activitiesInput: DG.InputBase;
    activitiesInputRoot: HTMLElement;
    similarityInput: DG.InputBase;

    constructor(editorSettings: DimReductionEditorOptions = {}) {
      super({...editorSettings, enableMCL: true});
      const numericalColumns = Array.from(this.tableInput.value!.columns.numericalNoDateTime);
      this.activitiesInput = ui.input.column('Activities', {table: this.tableInput.value!,
        value: numericalColumns.length ? numericalColumns[0] : null,
        filter: (col: DG.Column) => numericalColumns.includes(col)} as ColumnInputOptions);
      this.activitiesInputRoot = this.activitiesInput.root;
      this.similarityInput = ui.input.int('Similarity cutoff', {value: 80});
      this.methodsParams[MCLMethodName] = new MCLOptions() as any;
    }

    onTableInputChanged() {
      super.onTableInputChanged();
      if (this.activitiesInputRoot) {
        ui.empty(this.activitiesInputRoot);
        const numericalColumns = Array.from(this.tableInput.value!.columns.numerical);
        this.activitiesInput = ui.input.column('Activities', {table: this.tableInput.value!,
          value: numericalColumns.length ? numericalColumns[0] : null,
          filter: (col: DG.Column) => numericalColumns.includes(col)} as ColumnInputOptions);
        Array.from(this.activitiesInput.root.children).forEach((it) => this.activitiesInputRoot.append(it));
      }
    }

    public getEditor(): HTMLElement {
      return ui.div([
        this.tableInput,
        this.colInputRoot,
        this.preprocessingFunctionInputRoot,
        this.preprocessingFuncSettingsDiv,
        this.activitiesInputRoot,
        this.similarityMetricInputRoot,
        this.methodInput,
        this.methodSettingsDiv,
        this.similarityMetricInputRoot,
        this.similarityInput,
      ], {style: {minWidth: '320px'}, classes: 'ui-form dim-reduction-dialog-form'});
    }

    public getParams(): ActivityCliffsParams {
      return {
        ...super.getParams(),
        activities: this.activitiesInput.value!,
        similarityThreshold: this.similarityInput.value!
      };
    }

    public getInput() {
      return {...super.getInput(),
        activityCol: this.activitiesInput.value!.name, simThreshold: this.similarityInput.value};
    }

    public async applyInput(input: ReturnType<typeof this.getInput>) {
      super.applyInput(input);
      this.activitiesInput.value = this.tableInput.value!.col(input.activityCol);
      this.similarityInput.value = input.simThreshold;
    }
}
