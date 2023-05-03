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
        similarity: this.similarityInput.value!, methodName: this.methodInput.value!, options: this.algorithmOptions};
    }
  
    get paramsUI(): HTMLDivElement{
      return this.funcParamsDiv;
    }
  
    constructor(semtype: DG.SemType){
      super(semtype);
      this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, this.tableInput.value!.columns.byIndex(0));
      this.activitiesInputRoot = this.activitiesInput.root;
      this.similarityInput = ui.intInput('Similarity', 80);
      //@ts-ignore
      this.funcParamsDiv = ui.inputs([
        this.tableInput,
        this.molColInput,
        this.activitiesInput,
        this.similarityInput,
        this.methodInput,
        this.methodSettingsDiv
      ], {style: {minWidth: '320px'}});
    }

    onTableInputChanged(semtype: DG.SemType) {
        super.onTableInputChanged(semtype);
        ui.empty(this.activitiesInputRoot);
        this.activitiesInput = ui.columnInput('Activities', this.tableInput.value!, this.tableInput.value!.columns.byIndex(0));
        Array.from(this.activitiesInput.root.children).forEach((it) => this.activitiesInputRoot.append(it));
    }
  }
