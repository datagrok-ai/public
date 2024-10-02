/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import _ from 'lodash';

import {EventBus} from '../../model/event-bus';
import {PatternConfiguration, StrandType} from '../../model/types';
import {SubscriptionManager} from '../../model/subscription-manager';
import {STRAND, STRANDS, STRAND_LABEL, TERMINI, TERMINUS} from '../../model/const';
import {StringInput} from '../types';

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
    const controls = new StrandTerminalModificationControls(this.eventBus).create();

    $(editorBody).empty();
    $(editorBody).append(controls);
  }

  private resetToInitialState(): void {
    // this.eventBus.setLastLoadedPatternConfig(this.initialPatternConfig);
    this.eventBus.setPatternConfig(this.initialPatternConfig);
  }
}

class StrandTerminalModificationControls {
  constructor(
    private eventBus: EventBus
  ) { }

  create(): HTMLDivElement {
    const inputPanels = STRANDS.map((strand) => this.constructControlsPanel(strand));

    const container = ui.divH(inputPanels, {style: {gap: '24px'}});
    return container;
  }

  private constructControlsPanel(strand: StrandType): HTMLDivElement {
    if (!this.eventBus.isAntisenseStrandActive() && strand === STRAND.ANTISENSE)
      return ui.div([]);

    const modificationControls = this.createInputs(strand);

    const container = ui.block([
      ui.h1(`${STRAND_LABEL[strand]}`),
      modificationControls,
    ], {style: {paddingTop: '12px'}});

    return container;
  }

  private createInputs(strand: StrandType): HTMLElement {
    const termini = strand === STRAND.SENSE ? [...TERMINI].reverse() : TERMINI;
    const inputs = termini.map((terminus) => this.createInputForTerminus(strand, terminus)) as DG.InputBase[];
    return ui.form(inputs);
  }

  private createInputForTerminus(strand: StrandType, terminus: TERMINUS): StringInput {
    const initialValue = this.eventBus.getTerminalModifications()[strand][terminus];
    const input = ui.input.textArea(terminus, {value: initialValue});
    this.applyStylingToInput(input);

    input.onInput.subscribe(() => {
      const newValue = input.value;
      if (newValue === null)
        return;

      this.eventBus.updateTerminusModification(strand, terminus, newValue);
    });

    return input;
  }

  private applyStylingToInput(input: StringInput): void {
    const textarea = input.root.getElementsByTagName('textarea')[0];
    $(textarea).css('resize', 'none');
  }
}
