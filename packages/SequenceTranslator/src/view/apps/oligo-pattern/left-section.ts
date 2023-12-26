/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { SENSE_STRAND, ANTISENSE_STRAND, STRAND_LABEL, STRANDS, StrandType, OTHER_USERS } from '../../../model/pattern-app/const';

import {BooleanInput, StringInput, NumberInput} from './types';

import {EventBus} from '../../../model/pattern-app/event-bus';
import {PatternAppDataManager} from '../../../model/pattern-app/external-data-manager';
import {PatternConfigurationManager} from '../../../model/pattern-app/pattern-state-manager';
import $ from 'cash-dom';

export class PatternAppLeftSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
  ) {
    this.patternConfiguration = new PatternConfigurationManager(this.eventBus, this.dataManager);
  };
  private patternConfiguration: PatternConfigurationManager;

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
    private dataManager: PatternAppDataManager,
  ) { }

  getUiElements(): HTMLElement[] {
    return [
      ui.h1('Pattern'),
      this.toggleAntisenseStrandControl,
      this.strandLengthInputs[SENSE_STRAND].root,
      this.strandLengthInputs[ANTISENSE_STRAND].root,
      this.sequenceBaseInput.root,
      this.patternCommentInput.root,
      this.selectPattentInputBlock,
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
      `${STRAND_LABEL.ANTISENSE_STRAND} strand`, true, (visible: boolean) => this.eventBus.setAntisenseStrandVisibility(visible));
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

    this.eventBus.antisenseStrandActive$.subscribe((visible: boolean) => {
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

  private get selectPattentInputBlock(): HTMLDivElement {
    const patternChoiceControls = new PatternChoiceControls(
      this.eventBus,
      this.dataManager,
    );
    return patternChoiceControls.getControlsContainer();
  }
}

class PatternChoiceControls {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
  ) {
    this.eventBus.requestLoadPattern$.subscribe((value: string) => this.handlePatternChoice(value));

    const defaultUser = this.dataManager.getCurrentUserName();
    this.selectedUser = defaultUser;

    const defaultPattern = this.dataManager.getCurrentUserPatternNames()[0];
    this.selectedPattern = defaultPattern;

    this.patternChoiceContainer = ui.div([]);
    this.eventBus.patternListUpdate$.subscribe(() => this.updatePatternChoiceInputContainer()); 
  }

  private selectedUser: string;
  private selectedPattern: string;
  private patternChoiceContainer: HTMLDivElement;

  private handleUserChoice(userName: string) {
    this.selectedUser = userName;
    this.updatePatternChoiceInputContainer();
  }

  private handlePatternChoice(patternName: string) {
    this.selectedPattern = patternName;
    grok.shell.info(`Pattern ${patternName} selected`);
  }

  private isCurrentUserSelected(): boolean {
    return this.selectedUser !==  OTHER_USERS;
  }

  getControlsContainer(): HTMLDivElement {
    const patternInputs = this.getPatternInputs();
    this.patternChoiceContainer.append(patternInputs.root);
    return this.patternChoiceContainer;
  }

  private getPatternInputs(): StringInput {
    const userChoiceInput = this.userChoiceInput;
    const patternChoiceInput = this.getPatternChoiceInput();

    // todo: refactor this legacy solution
    patternChoiceInput.root.append(
      userChoiceInput.input,
      patternChoiceInput.input,
    );

    this.setPatternChoiceInputStyle(patternChoiceInput);

    const deletePatternButton = this.deletePatternButton;
    patternChoiceInput.addOptions(deletePatternButton);

    return patternChoiceInput;
  }

  private setPatternChoiceInputStyle(patternChoiceInput: StringInput): void {
    patternChoiceInput.setTooltip('Choose and apply pattern');
    patternChoiceInput.input.style.maxWidth = '120px';
    patternChoiceInput.input.style.marginLeft = '12px';
  }

  private get userChoiceInput(): StringInput {
    const currentUser = this.dataManager.getCurrentUserName();
    const possibleValues = [currentUser, OTHER_USERS];

    const userChoiceInput = ui.choiceInput(
      '', this.selectedUser, possibleValues,
      (userName: string) => this.handleUserChoice(userName)
    );
    this.setUserChoiceInputStyle(userChoiceInput);

    return userChoiceInput;
  }

  private setUserChoiceInputStyle(userChoiceInput: StringInput): void {
    userChoiceInput.setTooltip('Choose user to load pattern from');
    userChoiceInput.input.style.maxWidth = '142px';
  }

  private getPatternChoiceInput(): StringInput {
    const patternList = this.isCurrentUserSelected() ? this.dataManager.getCurrentUserPatternNames() : this.dataManager.getOtherUsersPatternNames();
    this.selectedPattern = patternList[0] || '';
    this.eventBus.requestPatternLoad(this.selectedPattern);
    const choiceInput = ui.choiceInput('Load pattern', this.selectedPattern, patternList, (value: string) => this.eventBus.requestPatternLoad(value));
    return choiceInput;
  }

  private get deletePatternButton(): HTMLDivElement {
    const button = ui.button(ui.iconFA('trash-alt'), () => this.eventBus.deletePattern(this.selectedPattern));

    return ui.div([ button ], 'ui-input-options');
  }

  private updatePatternChoiceInputContainer(): void {
    const patternInputs = this.getPatternInputs();
    $(this.patternChoiceContainer).empty();
    this.patternChoiceContainer.append(patternInputs.root);
  }
}
