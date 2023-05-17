import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IDimReductionParam, TSNEOptions, UMAPOptions } from '../reduce-dimensionality';
import { SequenceSpaceBaseFuncEditor } from './seq-space-base-editor';

export class ActivityCliffsFunctionEditor extends SequenceSpaceBaseFuncEditor {

    activitiesDataframe: DG.DataFrame | string;
    activitiesInput: DG.InputBase;
    activitiesInputRoot: HTMLElement;
    similarityInput: DG.InputBase;
    funcParamsDiv: HTMLDivElement;


    get funcParams(): any {
      return {table: this.tableInput.value!, molecules: this.molColInput.value!, activities: this.tableInput.value.col(this.activitiesDataframe ? this.activitiesInput.value!.name : '')!,
        similarity: this.similarityInput.value!, methodName: this.methodInput.value!, similarityMetric: this.similarityMetricInput.value!,
        options: this.algorithmOptions};
    }
  
    get paramsUI(): HTMLDivElement{
      return this.funcParamsDiv;
    }
  
    constructor(semtype: DG.SemType){
      super(semtype);
      const numericColumns = Array.from((this.tableInput.value as DG.DataFrame).columns.numerical);
      this.activitiesDataframe = numericColumns.length ? DG.DataFrame.fromColumns(numericColumns) : '';
      //@ts-ignore
      this.activitiesInput = ui.columnInput('Activities', this.activitiesDataframe, this.activitiesDataframe ? this.activitiesDataframe.columns.byIndex(0) : '');
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
        this.similarityMetricInput,
      ], {style: {minWidth: '320px'}});
    }

    onTableInputChanged(semtype: DG.SemType) {
        super.onTableInputChanged(semtype);
        ui.empty(this.activitiesInputRoot);
        const numericColumns = Array.from((this.tableInput.value as DG.DataFrame).columns.numerical);
        this.activitiesDataframe = numericColumns.length ? DG.DataFrame.fromColumns(numericColumns) : '';
        //@ts-ignore
        this.activitiesInput = ui.columnInput('Activities', this.activitiesDataframe, this.activitiesDataframe ? this.activitiesDataframe.columns.byIndex(0) : '');
        Array.from(this.activitiesInput.root.children).forEach((it) => this.activitiesInputRoot.append(it));
    }
  }
