/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import _ from 'lodash';

import {EventBus} from '../../../model/event-bus';
import {PatternConfiguration} from '../../../model/types';
import {HeaderControls} from './header-controls';
import {StrandControls} from './strand-controls';
import {SubscriptionManager} from '../../../model/subscription-manager';

export class StrandEditorDialog {
  private static isDialogOpen = false;
  private static instance: StrandEditorDialog;

  private initialPatternConfig: PatternConfiguration;
  private subscriptions = new SubscriptionManager();

  private constructor(
    private eventBus: EventBus
  ) { }

  static open(eventBus: EventBus): void {
    if (StrandEditorDialog.isDialogOpen)
      return;

    if (!StrandEditorDialog.instance)
      StrandEditorDialog.instance = new StrandEditorDialog(eventBus);

    StrandEditorDialog.instance.openDialog();
  }

  private openDialog(): void {
    this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
    StrandEditorDialog.isDialogOpen = true;
    this.createDialog().show();
  }

  private createDialog(): DG.Dialog {
    const editorBody = ui.divV([]);
    this.subscriptions.add(
      this.eventBus.strandsUpdated$.subscribe(() => this.onStrandsUpdated(editorBody))
    );

    const dialog = ui.dialog('Edit strands')
      .add(editorBody)
      .onOK(() => {})
      .onCancel(() => this.resetToInitialState());

    this.subscriptions.add(
      dialog.onClose.subscribe(() => {
        StrandEditorDialog.isDialogOpen = false;
        this.subscriptions.unsubscribeAll();
      })
    );

    return dialog;
  }

  private onStrandsUpdated(editorBody: HTMLDivElement) {
    const headerControls = new HeaderControls(
      this.eventBus, this.initialPatternConfig, this.subscriptions
    ).create();
    const strandControls = new StrandControls(this.eventBus, this.subscriptions).create();

    $(editorBody).empty();
    $(editorBody).append(headerControls, strandControls);
  }

  private resetToInitialState(): void {
    // this.eventBus.setLastLoadedPatternConfig(this.initialPatternConfig);
    this.eventBus.setPatternConfig(this.initialPatternConfig);
  }
}

