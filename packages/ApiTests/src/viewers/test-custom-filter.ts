import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

export class TestCustomFilter extends DG.Filter {
  get filterSummary(): string {
    return '';
  }

  constructor() {
    super();
  }

  applyState(state: any) {
    if (typeof state.column == 'string' || state.column instanceof String) {
      const columnName: string = state.column;
      delete state.column;
      state.columnName = columnName;
    }
    super.applyState(state);
  }

  // -- Filter --

  applyFilter(): void {
    const k = 11;
  }
}
