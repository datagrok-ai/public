/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { SENSE_STRAND, ANTISENSE_STRAND, STRAND_LABEL, STRANDS } from '../../../model/pattern-app/const';

import {BooleanInput, StringInput, NumberInput} from './types';

import {EventBus} from '../../../model/pattern-app/event-bus';
import {ExternalDataManager} from '../../../model/pattern-app/external-data-manager';
import {PatternConfigurationManager} from '../../../model/pattern-app/pattern-state-manager';

export class LeftSection {
  constructor(private eventBus: EventBus) {
    this.externalDataManager = ExternalDataManager.getInstance();
    this.patternConfiguration = new PatternConfigurationManager(this.eventBus);
  };
  private patternConfiguration: PatternConfigurationManager;
  private externalDataManager: ExternalDataManager;

  getLayout(): HTMLDivElement {
    const patternControlsManager = new PatternControlsManager(
      this.eventBus,
      this.patternConfiguration,
      this.externalDataManager
    );
    const patternConstrolsBlock = patternControlsManager.getUiElements();
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
  constructor(
    private eventBus: EventBus,
    private patternConfiguration: PatternConfigurationManager,
    private externalDataManager: ExternalDataManager,
  ) { }

  getUiElements(): HTMLElement[] {
    return [
      ui.h1('Pattern'),
      this.toggleAntisenseStrandControl,
      this.strandLengthInputs[SENSE_STRAND].root,
      this.strandLengthInputs[ANTISENSE_STRAND].root,
      this.sequenceBaseInput.root,
      this.patternCommentInput.root,
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
    const strandLengthInputs = Object.fromEntries(
      STRANDS.map(
        (strand) => {
          const sequenceLength = this.patternConfiguration.getBases(strand).length;
          const input = ui.intInput(`${STRAND_LABEL[strand]} length`, sequenceLength);
          input.setTooltip(`Length of ${STRAND_LABEL[strand].toLowerCase()}, including overhangs`);
          return [strand, input];
        }
      )
    );

    this.eventBus.antisenseStrandVisible$.subscribe((visible: boolean) => {
      $(strandLengthInputs[ANTISENSE_STRAND].root).toggle(visible);
    })

    return strandLengthInputs;
  }

  private get sequenceBaseInput(): StringInput {
    const nucleotideBaseChoices = this.externalDataManager.fetchNucleotideBases();
    const defaultNucleotideBase = nucleotideBaseChoices[0];

    const sequenceBaseInput = ui.choiceInput('Sequence basis', defaultNucleotideBase, nucleotideBaseChoices, (value: string) => {
    });
    sequenceBaseInput.setTooltip('Nucleotide base to use for the sequence');
    return sequenceBaseInput;
  }

  private get patternCommentInput(): StringInput {
    const patternCommentInput = ui.textInput('Comment', '', (value: string) => this.eventBus.changeComment(value));
    return patternCommentInput;
  }
}
