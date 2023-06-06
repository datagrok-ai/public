import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IDimReductionParam, TSNEOptions, UMAPOptions } from '../reduce-dimensionality';
import { SequenceSpaceBaseFuncEditor } from './seq-space-base-editor';

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
      this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, this.tableInput.value!.columns.byIndex(0), null, {'predicate': (col: DG.Column) => col.type === DG.TYPE.INT});
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
        this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, this.tableInput.value!.columns.byIndex(0), null, {'predicate': (col: DG.Column) => col.type === DG.TYPE.INT});
        Array.from(this.activitiesInput.root.children).forEach((it) => this.activitiesInputRoot.append(it));
    }
  }