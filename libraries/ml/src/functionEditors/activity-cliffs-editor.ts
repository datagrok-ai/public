import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IDimReductionParam, TSNEOptions, UMAPOptions } from '../reduce-dimensionality';
import { SequenceSpaceBaseFuncEditor } from './seq-space-base-editor';
import { ColumnInputOptions } from '@datagrok-libraries/utils/src/type-declarations';

export class ActivityCliffsFunctionEditor extends SequenceSpaceBaseFuncEditor {

    activitiesInput: DG.InputBase;
    activitiesInputRoot: HTMLElement;
    similarityInput: DG.InputBase;
    funcParamsDiv: HTMLDivElement;


    get funcParams(): any {
      return {table: this.tableInput.value!, molecules: this.molColInput.value!, activities: this.activitiesInput.value!,
        similarity: this.similarityInput.value!, methodName: this.methodInput.value!, similarityMetric: this.similarityMetricInput.value!,
        options: this.algorithmOptions};
    }
  
    get paramsUI(): HTMLDivElement{
      return this.funcParamsDiv;
    }
  
    constructor(semtype: DG.SemType){
      super(semtype);
      const numericalColumns = this.tableInput.value!.columns.numerical;
      //TODO: remove when the new version of datagrok-api is available
      //@ts-ignore
      this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, DG.Utils.firstOrNull(numericalColumns), null, {filter: (col: DG.Column) => Array.from(numericalColumns).includes(col)} as ColumnInputOptions);
      this.activitiesInputRoot = this.activitiesInput.root;
      this.similarityInput = ui.intInput('Similarity cutoff', 80);
      ui.tooltip.bind(this.similarityInput.root, `Pairs of similar (cutoff is used) molecules with high difference in activity are considered 'cliffs'`)
      //@ts-ignore
      this.funcParamsDiv = ui.inputs([
        this.tableInput,
        this.molColInput,
        this.activitiesInput,
        this.similarityInput,
        this.methodInput,
        this.methodSettingsDiv,
        this.similarityMetricInputRoot,
      ], {style: {minWidth: '320px'}});
    }

    onTableInputChanged(semtype: DG.SemType) {
        super.onTableInputChanged(semtype);
        ui.empty(this.activitiesInputRoot);
        const numericalColumns = this.tableInput.value!.columns.numerical;
        //TODO: remove when the new version of datagrok-api is available
        //@ts-ignore
        this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, DG.Utils.firstOrNull(numericalColumns), null, {filter: (col: DG.Column) => Array.from(numericalColumns).includes(col)} as ColumnInputOptions);
        Array.from(this.activitiesInput.root.children).forEach((it) => this.activitiesInputRoot.append(it));
    }
  }