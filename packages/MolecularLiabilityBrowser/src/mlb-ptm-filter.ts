import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

export class MlbPtmFilter extends DG.Filter {
  get filterSummary(): string {
    return '';
  }

  constructor() {
    super();
  }

  applyState(state: any) {
    super.applyState(state);

    this.render();
  }

  // -- View --

  render(): void {
    $(this.root).empty();

    this.root.appendChild(ui.div('Hallo'));
  }

  // -- Filter --

  applyFilter(): void {
  }
}
