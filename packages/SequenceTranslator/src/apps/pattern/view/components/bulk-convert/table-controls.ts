import * as ui from 'datagrok-api/ui';
import '../../style.css';

import {EventBus} from '../../../model/event-bus';
import {ColumnInputManager} from './column-input';
import {TableInputManager} from './table-input';
import {STRAND} from '../../../model/const';

export class TableControlsManager {
  private tableInputManager: TableInputManager;
  private columnInputManager: ColumnInputManager;

  constructor(private eventBus: EventBus) {
    this.tableInputManager = new TableInputManager(eventBus);
    this.columnInputManager = new ColumnInputManager(eventBus);
  }

  createControls(): HTMLElement[] {
    const title = ui.h1('Bulk convert');
    const tableInput = this.tableInputManager.getTableInputContainer();
    const columnControls = this.columnInputManager.getColumnControlsContainer();

    const convertSequenceButton = ui.bigButton('Convert', () => this.processConvertButtonClick());

    return [
      title,
      tableInput,
      columnControls,
      ui.buttonsInput([
        convertSequenceButton,
      ]),
    ];
  }

  private processConvertButtonClick(): void {
    // grok.shell.info(`Selected table: ${this.eventBus.getTableSelection()?.name}, selected AS strand column: ${this.eventBus.getSelectedColumn(STRAND.ANTISENSE)}`);
  }
}
