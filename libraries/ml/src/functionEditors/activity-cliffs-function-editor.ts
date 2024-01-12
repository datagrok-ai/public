import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ColumnInputOptions} from '@datagrok-libraries/utils/src/type-declarations';
import {DimReductionBaseEditor, DimReductionEditorOptions, DimReductionParams} from './dimensionality-reduction-editor';

export type ActivityCliffsParams = DimReductionParams & {
    activities: DG.Column;
    similarityThreshold: number;
}

export class ActivityCliffsEditor extends DimReductionBaseEditor {
    activitiesInput: DG.InputBase;
    activitiesInputRoot: HTMLElement;
    similarityInput: DG.InputBase;

    constructor(editorSettings: DimReductionEditorOptions = {}) {
      super(editorSettings);
      const numericalColumns = this.tableInput.value!.columns.numerical;
      this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!,
        DG.Utils.firstOrNull(numericalColumns), null,
        {filter: (col: DG.Column) => Array.from(numericalColumns).includes(col)} as ColumnInputOptions);
      this.activitiesInputRoot = this.activitiesInput.root;
      this.similarityInput = ui.intInput('Similarity cutoff', 80);
    }

    onTableInputChanged() {
      super.onTableInputChanged();
      if (this.activitiesInputRoot) {
        ui.empty(this.activitiesInputRoot);
        const numericalColumns = this.tableInput.value!.columns.numerical;
        this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!,
          DG.Utils.firstOrNull(numericalColumns), null,
        {filter: (col: DG.Column) => Array.from(numericalColumns).includes(col)} as ColumnInputOptions);
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
      ], {style: {minWidth: '320px'}, classes: 'ui-form'});
    }

    public getParams(): ActivityCliffsParams {
      return {
        ...super.getParams(),
        activities: this.activitiesInput.value!,
        similarityThreshold: this.similarityInput.value!
      };
    }
}
