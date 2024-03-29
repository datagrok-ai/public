/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import '../style.css';

import {MAX_SEQUENCE_LENGTH, STRAND, STRANDS, STRAND_LABEL} from '../../model/const';
import {PatternDefaultsProvider} from '../../model/defaults-provider';
import {EventBus} from '../../model/event-bus';
import {StrandType} from '../../model/types';
import {NumberInput, StringInput} from '../types';
import {isDialogOpen, StrandEditorDialog} from './strand-editor/dialog';

export class PatternEditControlsManager {
  constructor(
    private eventBus: EventBus,
    private defaultState: PatternDefaultsProvider
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

    return [
      ui.h1('Edit'),
      antisenseStrandToggle,
      senseStrandLengthInput,
      antisenseStrandLengthInput,
      sequenceBaseInput,
      patternNameInputBlock,
      patternCommentInput,
      ui.buttonsInput([
        editPatternButton,
      ]),
    ];
  }

  private createEditPatternButton(): HTMLButtonElement {
    const editPatternButton = ui.button(
      'Edit strands',
      () => {
        if (!isDialogOpen)
          new StrandEditorDialog(this.eventBus).open();
      });

    ui.tooltip.bind(editPatternButton, 'Edit modifications and PTOs per strand');

    return editPatternButton;
  }


  private createAntisenseStrandToggle(): HTMLElement {
    const toggleAntisenseStrand = ui.switchInput(
      `${STRAND_LABEL[STRAND.ANTISENSE]} strand`,
      true
    );

    toggleAntisenseStrand.onInput(
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

      const input = ui.intInput(
        `${STRAND_LABEL[strand]} length`,
        sequenceLength
      );
      input.onInput(() => updateStrandLengthInputs(strand, input));

      this.eventBus.patternStateChanged$.subscribe(() => {
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
    const availableNucleoBases = this.defaultState.fetchAvailableNucleotideBases()
      .sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
    const defaultNucleotideBase = this.defaultState.fetchDefaultNucleobase();

    const sequenceBaseInput = ui.choiceInput('Sequence basis', defaultNucleotideBase, availableNucleoBases);

    sequenceBaseInput.onInput(() => this.eventBus.replaceSequenceBase(sequenceBaseInput.value!));

    this.eventBus.patternStateChanged$.subscribe(() => {
      sequenceBaseInput.value = this.eventBus.getSequenceBase();
    });

    sequenceBaseInput.setTooltip('Most frequent nucleobase in the strands');
    return sequenceBaseInput;
  }

  private createPatternCommentInput(): StringInput {
    const patternCommentInput = ui.textInput(
      'Comment',
      this.eventBus.getComment()
    );

    $(patternCommentInput.root).addClass('st-pattern-text-input');

    patternCommentInput.onInput(
      () => this.eventBus.updateComment(patternCommentInput.value!)
    );

    this.eventBus.patternLoaded$.subscribe(() => {
      patternCommentInput.value = this.eventBus.getComment();
    });

    return patternCommentInput;
  }

  private createPatternNameInputBlock(): HTMLElement {
    const patternNameInput = ui.textInput(
      'Pattern name',
      this.eventBus.getPatternName()
    );

    $(patternNameInput.root).addClass('st-pattern-text-input');

    patternNameInput.onInput(
      () => this.eventBus.updatePatternName(patternNameInput.value)
    );
    this.eventBus.patternLoaded$.subscribe(() => {
      patternNameInput.value = this.eventBus.getPatternName();
    });

    patternNameInput.setTooltip('Name under which pattern will be saved');

    return patternNameInput.root;
  }
}

