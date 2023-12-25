/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { SENSE_STRAND, ANTISENSE_STRAND, STRAND_LABEL, STRANDS, StrandType, OTHER_USERS } from '../../../model/pattern-app/const';

import {BooleanInput, StringInput, NumberInput} from './types';

import {EventBus} from '../../../model/pattern-app/event-bus';
import {AppDataManager} from '../../../model/pattern-app/external-data-manager';
import {PatternConfigurationManager} from '../../../model/pattern-app/pattern-state-manager';

export class LeftSection {
  constructor(private eventBus: EventBus) {
    this.dataManager = AppDataManager.getInstance(this.eventBus);
    this.patternConfiguration = new PatternConfigurationManager(this.eventBus);
  };
  private patternConfiguration: PatternConfigurationManager;
  private dataManager: AppDataManager;

  getLayout(): HTMLDivElement {
    const patternControlsManager = new PatternControlsManager(
      this.eventBus,
      this.patternConfiguration,
      this.dataManager
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
    private dataManager: AppDataManager,
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
    const createStrandLengthInput = (strand: StrandType) => {
      const sequenceLength = this.patternConfiguration.getBases(strand).length;
      const input = ui.intInput(`${STRAND_LABEL[strand]} length`, sequenceLength);
      input.setTooltip(`Length of ${STRAND_LABEL[strand].toLowerCase()}, including overhangs`);
      return [strand, input];
    }

    const strandLengthInputs = Object.fromEntries(
      STRANDS.map((strand) => createStrandLengthInput(strand))
    );

    this.eventBus.antisenseStrandVisible$.subscribe((visible: boolean) => {
      $(strandLengthInputs[ANTISENSE_STRAND].root).toggle(visible);
    })

    return strandLengthInputs;
  }

  private get sequenceBaseInput(): StringInput {
    const nucleotideBaseChoices = this.dataManager.fetchNucleotideBases();
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

class PatternChoiceControls {
  constructor(
    private eventBus: EventBus,
    // private patternConfiguration: PatternConfigurationManager,
    private dataManager: AppDataManager,
  ) {
    this.eventBus.userChoice$.subscribe((value: string) => this.handleUserChoice(value));
    this.eventBus.loadPattern$.subscribe((value: string) => this.handlePatternChoice(value));

    const defaultUser = this.dataManager.getCurrentUserName();
    this.selectedUser = defaultUser;

    const defaultPattern = this.dataManager.getCurrentUserPatterns()[0];
    this.selectedPattern = defaultPattern;
  }

  private selectedUser: string;
  private selectedPattern: string;

  private handleUserChoice(value: string) {
    this.selectedUser = value;
    grok.shell.info(`User ${value} selected`);
  }

  private handlePatternChoice(value: string) {
    this.selectedPattern = value;
    grok.shell.info(`Pattern ${value} selected`);
  }

  private isCurrentUserSelected(): boolean {
    return this.selectedUser !==  OTHER_USERS;
  }

  private get userChoiceInput(): StringInput {
    const currentUser = this.dataManager.getCurrentUserName();

    const values = [currentUser, OTHER_USERS];
    const userChoiceInput = ui.choiceInput(
      '', currentUser, values,
      (value: string) => this.eventBus.chooseUser(value)
    );
    userChoiceInput.setTooltip('Choose user to load pattern from');

    return userChoiceInput;
  }

  private getPatternChoiceInput(): StringInput {
    const patternList = this.isCurrentUserSelected() ? this.dataManager.getCurrentUserPatterns() : this.dataManager.getOtherUsersPatterns();
    const choiceInput = ui.choiceInput('Load pattern', this.selectedPattern, patternList, (value: string) => this.eventBus.loadPattern(value));
    return choiceInput;
  }

  private get deletePatternButton(): HTMLDivElement {
    const button = ui.button(ui.iconFA('trash-alt'), () => this.eventBus.deletePattern(this.selectedPattern));

    return ui.div([ button, ], 'ui-input-options');
  }
}
