import * as ui from 'datagrok-api/ui';
import '../../style.css';

import {EventBus} from '../../../model/event-bus';
import {ColumnInputManager} from './column-input';
import {TableInputManager} from './table-input';
import {bulkTranslate} from '../../../model/translator';

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

    const convertButton = ui.bigButton('Convert', () => this.processConvertButtonClick());

    return [
      title,
      tableInput,
      columnControls,
      ui.buttonsInput([
        convertButton,
      ]),
    ];
  }

  private processConvertButtonClick(): void {
    bulkTranslate(this.eventBus);
  }
}
