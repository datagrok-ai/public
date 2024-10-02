/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import '../style.css';

import {MAX_SEQUENCE_LENGTH, STRAND, STRANDS, STRAND_LABEL} from '../../model/const';
import {EventBus} from '../../model/event-bus';
import {StrandType} from '../../model/types';
import {NumberInput, StringInput} from '../types';
import {StrandEditorDialog} from './strand-editor/dialog';
import {DataManager} from '../../model/data-manager';
import {TerminalModificationEditorDialog} from './terminal-modification-editor';

export class PatternEditControlsManager {
  constructor(
    private eventBus: EventBus,
    private dataManager: DataManager
  ) { }

  createControls(): HTMLElement[] {
    const antisenseStrandToggle = this.createAntisenseStrandToggle();
    const strandLengthInputs = this.createStrandLengthInputs();

    const senseStrandLengthInput = strandLengthInputs[STRAND.SENSE].root;
    const antisenseStrandLengthInput = strandLengthInputs[STRAND.ANTISENSE].root;

    const sequenceBaseInput = this.createSequenceBaseInput().root;
    const patternCommentInput = this.createPatternCommentInput().root;
    const patternNameInputBlock = this.createPatternNameInputBlock();

    const editPatternButton = this.createEditPatternButton();
    const editTerminalModificationsButton = this.createEditTerminalModificationsButton();

    return [
      ui.h1('Edit'),
      antisenseStrandToggle,
      senseStrandLengthInput,
      antisenseStrandLengthInput,
      sequenceBaseInput,
      patternNameInputBlock,
      patternCommentInput,
      ui.buttonsInput([
        editTerminalModificationsButton,
        editPatternButton,
      ]),
    ];
  }

  private createEditPatternButton(): HTMLButtonElement {
    const editPatternButton = ui.button(
      'Edit strands',
      () => StrandEditorDialog.open(this.eventBus, this.dataManager)
    );

    ui.tooltip.bind(editPatternButton, 'Edit strand modifications and PTOs');

    return editPatternButton;
  }

  private createEditTerminalModificationsButton(): HTMLButtonElement {
    const editTerminalModificationsButton = ui.button(
      'Edit terminals',
      () => TerminalModificationEditorDialog.open(this.eventBus)
    );

    ui.tooltip.bind(editTerminalModificationsButton, 'Edit terminal modifications');
    $(editTerminalModificationsButton).css('margin-right', '20px');

    return editTerminalModificationsButton;
  }


  private createAntisenseStrandToggle(): HTMLElement {
    const toggleAntisenseStrand = ui.input.toggle(
      `${STRAND_LABEL[STRAND.ANTISENSE]} strand`,
      {value: this.eventBus.isAntisenseStrandActive()}
    );

    toggleAntisenseStrand.onInput.subscribe(
      () => this.eventBus.toggleAntisenseStrand(toggleAntisenseStrand.value)
    );

    this.eventBus.patternLoaded$.subscribe(() => {
      toggleAntisenseStrand.value = this.eventBus.isAntisenseStrandActive();
    });

    toggleAntisenseStrand.setTooltip('Toggle antisense strand');

    return toggleAntisenseStrand.root;
  }

  private createStrandLengthInputs(): Record<string, NumberInput> {
    const createStrandLengthInput = (strand: StrandType) => {
      const sequenceLength = this.eventBus.getNucleotideSequences()[strand].length;

      const input = ui.input.int(`${STRAND_LABEL[strand]} length`, {value: sequenceLength});
      input.onInput.subscribe(() => updateStrandLengthInputs(strand, input));

      this.eventBus.nucleotideSequencesChanged$.subscribe(() => {
        input.value = this.eventBus.getNucleotideSequences()[strand].length;
      });

      input.setTooltip(`Number of nucleotides in ${strand}, including overhangs`);
      return [strand, input];
    };

    const updateStrandLengthInputs = (strand: StrandType, input: DG.InputBase<number | null>) => {
      const length = input.value;
      if (length === null) return;

      if (length <= 0) {
        grok.shell.warning(`Sequence length must be greater than 0`);
        input.value = 1;
      }
      if (length > MAX_SEQUENCE_LENGTH) {
        grok.shell.warning(`Sequence length must be less than ${MAX_SEQUENCE_LENGTH + 1}`);
        input.value = MAX_SEQUENCE_LENGTH;
      }

      this.eventBus.updateStrandLength(strand, input.value!);
    };

    const strandLengthInputs = Object.fromEntries(
      STRANDS.map((strand) => createStrandLengthInput(strand))
    );

    this.eventBus.antisenseStrandToggled$.subscribe((active: boolean) => {
      $(strandLengthInputs[STRAND.ANTISENSE].root).toggle(active);
    });

    return strandLengthInputs;
  }

  private createSequenceBaseInput(): StringInput {
    const availableNucleoBases = this.dataManager.fetchAvailableNucleotideBases()
      .sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
    const defaultNucleotideBase = this.dataManager.fetchDefaultNucleobase();

    const sequenceBaseInput = ui.input.choice('Sequence basis', {value: defaultNucleotideBase, items: availableNucleoBases});

    sequenceBaseInput.onInput.subscribe(() => this.eventBus.replaceSequenceBase(sequenceBaseInput.value!));

    this.eventBus.nucleotideSequencesChanged$.subscribe(() => {
      sequenceBaseInput.value = this.eventBus.getSequenceBase();
    });

    sequenceBaseInput.setTooltip('Most frequent nucleobase in the strands');
    return sequenceBaseInput;
  }

  private createPatternCommentInput(): StringInput {
    const patternCommentInput = ui.input.textArea('Comment', {value: this.eventBus.getComment()});

    $(patternCommentInput.root).addClass('st-pattern-text-input');

    patternCommentInput.onInput.subscribe(
      () => this.eventBus.updateComment(patternCommentInput.value!)
    );

    this.eventBus.patternLoaded$.subscribe(() => {
      patternCommentInput.value = this.eventBus.getComment();
    });

    return patternCommentInput;
  }

  private createPatternNameInputBlock(): HTMLElement {
    const patternNameInput = ui.input.textArea('Pattern name', {value: this.eventBus.getPatternName()});

    $(patternNameInput.root).addClass('st-pattern-text-input');

    patternNameInput.onInput.subscribe(
      () => this.eventBus.updatePatternName(patternNameInput.value)
    );
    this.eventBus.patternLoaded$.subscribe(() => {
      patternNameInput.value = this.eventBus.getPatternName();
    });

    patternNameInput.setTooltip('Name under which pattern will be saved');

    return patternNameInput.root;
  }
}

