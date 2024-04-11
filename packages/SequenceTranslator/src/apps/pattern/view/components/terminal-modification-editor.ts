/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import _ from 'lodash';

import {EventBus} from '../../model/event-bus';
import {PatternConfiguration} from '../../model/types';
import {SubscriptionManager} from '../../model/subscription-manager';

export class TerminalModificationEditorDialog {
  private static isDialogOpen = false;
  private static instance: TerminalModificationEditorDialog;

  private initialPatternConfig: PatternConfiguration;
  private subscriptions = new SubscriptionManager();

  private constructor(
    private eventBus: EventBus
  ) { }

  static open(eventBus: EventBus): void {
    if (TerminalModificationEditorDialog.isDialogOpen)
      return;

    if (!TerminalModificationEditorDialog.instance)
      TerminalModificationEditorDialog.instance = new TerminalModificationEditorDialog(eventBus);

    TerminalModificationEditorDialog.instance.openDialog();
  }

  private openDialog(): void {
    this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
    TerminalModificationEditorDialog.isDialogOpen = true;
    this.createDialog().show();
  }

  private createDialog(): DG.Dialog {
    const editorBody = ui.divV([]);
    this.subscriptions.add(
      this.eventBus.strandsUpdated$.subscribe(() => this.onStrandsUpdated(editorBody))
    );

    const dialog = ui.dialog('Edit terminal modifications')
      .add(editorBody)
      .onOK(() => {})
      .onCancel(() => this.resetToInitialState());

    this.subscriptions.add(
      dialog.onClose.subscribe(() => {
        TerminalModificationEditorDialog.isDialogOpen = false;
        this.subscriptions.unsubscribeAll();
      })
    );

    return dialog;
  }

  private onStrandsUpdated(editorBody: HTMLDivElement) {
    this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
    const controls = ui.divText('Controls');

    $(editorBody).empty();
    $(editorBody).append(controls);
  }

  private resetToInitialState(): void {
    this.eventBus.setPatternConfig(this.initialPatternConfig);
  }
}

