/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULT_PHOSPHOROTHIOATE, DEFAULT_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, SENSE_STRAND, ANTISENSE_STRAND, STRAND_LABEL, STRANDS, TERMINAL_KEYS, TERMINAL, THREE_PRIME_END, FIVE_PRIME_END, PATTERN_FIELD as PATTERN_KEY, StrandType, TerminalType
} from '../../../model/pattern-app/const';

import {PatternConfiguration, BooleanInput, StringInput, NumberInput} from './types';
import {applyToAllStrands} from './utils';

import {EventBus} from '../../../model/pattern-app/event-bus';

export class LeftSection {
  private eventBus = new EventBus();

  getLayout(): HTMLDivElement {
    const patternConstrolsBlock = new PatternControlsManager(this.eventBus).getUiElements();
    const layout = ui.box(
      ui.div([
          ...patternConstrolsBlock
        ],
        'ui-form'
      ),
      {style: {maxWidth: '450px'}}
    );
    return layout;
  }
}

export class PatternControlsManager {
  constructor(private eventBus: EventBus) { }

  getUiElements(): HTMLElement[] {
    return [
      ui.h1('Pattern'),
      this.toggleAntisenseStrandControl,
      this.strandLengthInputs[SENSE_STRAND].root,
      this.strandLengthInputs[ANTISENSE_STRAND].root,
      // sequenceBase.root,
      // patternCommentInput.root,
      // loadPatternDiv,
      // patternNameInput.root,
      // ui.h1('Convert'),
      // tableInput.root,
      // strandColumnInput[SENSE_STRAND],
      // strandColumnInput[ANTISENSE_STRAND],
      // idColumnSelector.root,
      // ui.buttonsInput([convertSequenceButton]),
    ]
  }

  private get toggleAntisenseStrandControl(): HTMLElement {
    const toggleAntisenseStrand = ui.switchInput(
      `${STRAND_LABEL.ANTISENSE_STRAND} strand`, true, (value: boolean) => this.eventBus.toggleAntisenseStrand(value));
    toggleAntisenseStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    return toggleAntisenseStrand.root;
  }

  private get strandLengthInputs(): Record<string, NumberInput> {
    const strandLengthInputs = applyToAllStrands(
      (strand) => {
        const input = ui.intInput(`${STRAND_LABEL[strand]} length`, DEFAULT_SEQUENCE_LENGTH);
        input.setTooltip(`Length of ${STRAND_LABEL[strand].toLowerCase()}, including overhangs`);
        return [strand, input];
      }
    );

    this.eventBus.antisenseStrandVisible$.subscribe((createAsCriterion: boolean) => {
      $(strandLengthInputs[ANTISENSE_STRAND].root).toggle(createAsCriterion);
    })

    return strandLengthInputs;
  }
}
