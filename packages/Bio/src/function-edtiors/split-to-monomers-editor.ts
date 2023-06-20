import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class SplitToMonomersFunctionEditor {
  tableInput: DG.InputBase;
  seqColInput: DG.InputBase;

  funcParamsDiv: HTMLDivElement;

  get funcParams(): {} {
    return {
      table: this.tableInput.value!,
      sequence: this.seqColInput.value!,
    };
  }

  get paramsUI(): HTMLDivElement {
    return this.funcParamsDiv;
  }

  constructor() {
    this.tableInput = ui.tableInput('Table', grok.shell.tv.dataFrame, undefined, () => {
      this.onTableInputChanged();
    });
    //TODO: remove when the new version of datagrok-api is available
    const seqColValue = this.tableInput.value!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    const seqColOptions = {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE};
    //@ts-ignore
    this.seqColInput = ui.columnInput('Sequence', this.tableInput.value!, seqColValue, null, seqColOptions);

    this.funcParamsDiv = ui.inputs([
      this.tableInput,
      this.seqColInput,
    ], {style: {minWidth: '320px'}});
  }

  onTableInputChanged(): void {
    this.seqColInput = ui.columnInput('Sequence', this.tableInput.value!,
      this.tableInput.value!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE));
  }
}
