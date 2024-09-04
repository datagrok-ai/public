import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export class InputColumnBase extends DG.InputBase<DG.Column | null> {

  private readonly options?: ui.input.IColumnInputInitOptions<DG.Column>;

  constructor(name: string, options?: ui.input.IColumnInputInitOptions<DG.Column>) {
    const base = ui.input.column(name, options);
    super(base.dart, options?.onValueChanged);
    this.options = options;
  }

  setColumnInputTable(table?: DG.DataFrame): void {
    ui.input.setColumnInputTable(this, table!, this.options?.filter);
  }
}

declare module 'datagrok-api/ui' {
  namespace input {
    function column2(name: string, options?: ui.input.IColumnInputInitOptions<DG.Column>): InputColumnBase;
  }
}

ui.input.column2 = function(name: string, options?: ui.input.IColumnInputInitOptions<DG.Column>): InputColumnBase {
  const res = new InputColumnBase(name, options);
  return res;
};
