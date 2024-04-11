/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import _ from 'lodash';

import {EventBus} from '../../../model/event-bus';
import {PatternConfiguration} from '../../../model/types';
import {HeaderControls} from './header-controls';
import {StrandControls} from './strand-controls';

export let isDialogOpen = false;

export class StrandEditorDialog {
  private initialPatternConfig: PatternConfiguration;

  constructor(
    private eventBus: EventBus
  ) {
    this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
  }

  open(): void {
    isDialogOpen = true;
    this.createDialog().show();
  }

  private createDialog(): DG.Dialog {
    const editorBody = ui.divV([]);
    this.eventBus.strandsUpdated$.subscribe(() => {
      this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
      const header = new HeaderControls(this.eventBus, this.initialPatternConfig).getPhosphorothioateLinkageControls();
      const controls = new StrandControls(this.eventBus).create();

      $(editorBody).empty();
      $(editorBody).append(header, controls);
    });

    const dialog = ui.dialog('Edit strands')
      .add(editorBody)
      .onOK(() => {})
      .onCancel(() => this.resetToInitialState());

    dialog.onClose.subscribe(() => isDialogOpen = false);

    return dialog;
  }

  private resetToInitialState(): void {
    this.eventBus.setPatternConfig(this.initialPatternConfig);
  }
}

