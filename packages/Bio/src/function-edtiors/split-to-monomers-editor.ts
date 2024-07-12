import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class SplitToMonomersFunctionEditor {
  tableInput: DG.InputBase;
  seqColInput: DG.InputBase;

  funcParamsDiv: HTMLElement;

  get funcParams(): {} {
    return {
      table: this.tableInput.value!,
      sequence: this.seqColInput.value!,
    };
  }

  get paramsUI(): HTMLElement {
    return this.funcParamsDiv;
  }

  constructor() {
    this.tableInput = ui.input.table('Table', {value: grok.shell.tv.dataFrame, onValueChanged: () => {
      this.onTableInputChanged();
    }});
    //TODO: remove when the new version of datagrok-api is available
    const seqColValue = this.tableInput.value!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    this.seqColInput = ui.input.column('Sequence', {table: this.tableInput.value!, value: seqColValue,
      filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE});

    this.funcParamsDiv = ui.inputs([
      this.tableInput,
      this.seqColInput,
    ], {style: {minWidth: '320px'}});
  }

  onTableInputChanged(): void {
    this.seqColInput = ui.input.column('Sequence', {table: this.tableInput.value!,
      value: this.tableInput.value!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)});
  }
}
