/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { SENSE_STRAND, ANTISENSE_STRAND, STRAND_LABEL, STRANDS, StrandType, OTHER_USERS } from '../model/const';

import {StringInput, NumberInput} from './types';

import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {PatternConfigurationManager} from '../model/pattern-state-manager';
import $ from 'cash-dom';

export class PatternAppLeftSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
    private patternConfiguration: PatternConfigurationManager,
  ) {
  };

  getLayout(): HTMLDivElement {
    const patternControlsManager = new PatternControlsManager(
      this.eventBus,
      this.patternConfiguration,
      this.dataManager
    );
    const tableControlsManager = new TableControlsManager(
      this.eventBus,
      this.patternConfiguration,
      this.dataManager
    );

    const patternConstrolsBlock = patternControlsManager.createUIComponents();
    const tableControlsBlock = tableControlsManager.createUIComponents();

    const layout = ui.box(
      ui.div([
          ...patternConstrolsBlock,
          ...tableControlsBlock
        ],
        'ui-form'
      ),
      {style: {maxWidth: '450px'}}
    );
    return layout;
  }
}

class PatternControlsManager {
  constructor(
    private eventBus: EventBus,
    private patternConfiguration: PatternConfigurationManager,
    private dataManager: PatternAppDataManager,
  ) { }

  createUIComponents(): HTMLElement[] {
    const title = ui.h1('Pattern');

    const antisenseStrandToggle = this.createAntisenseStrandToggle();
    const strandLengthInputs = this.createStrandLengthInputs();

    const senseStrandLengthInput = strandLengthInputs[SENSE_STRAND].root;
    const antisenseStrandLengthInput = strandLengthInputs[ANTISENSE_STRAND].root;

    const sequenceBaseInput = this.createSequenceBaseInput().root;
    const patternCommentInput = this.createPatternCommentInput().root;
    const patternSelectionBlock = this.createPatternSelectionBlock();
    const patternNameInputBlock = this.createPatternNameInputBlock();

    return [
      title,
      antisenseStrandToggle,
      senseStrandLengthInput,
      antisenseStrandLengthInput,
      sequenceBaseInput,
      patternCommentInput,
      patternSelectionBlock,
      patternNameInputBlock,
    ];
  }

  private createAntisenseStrandToggle(): HTMLElement {
    const toggleAntisenseStrand = ui.switchInput(
      `${STRAND_LABEL.ANTISENSE_STRAND} strand`, true, (isActive: boolean) => this.eventBus.toggleAntisenseStrand(isActive));
    toggleAntisenseStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    return toggleAntisenseStrand.root;
  }

  private createStrandLengthInputs(): Record<string, NumberInput> {
    const createStrandLengthInput = (strand: StrandType) => {
      const sequenceLength = this.patternConfiguration.getBases(strand).length;
      const input = ui.intInput(`${STRAND_LABEL[strand]} length`, sequenceLength);
      input.setTooltip(`Length of ${STRAND_LABEL[strand].toLowerCase()}, including overhangs`);
      return [strand, input];
    }

    const strandLengthInputs = Object.fromEntries(
      STRANDS.map((strand) => createStrandLengthInput(strand))
    );

    this.eventBus.isAntisenseStrandActive$.subscribe((active: boolean) => {
      $(strandLengthInputs[ANTISENSE_STRAND].root).toggle(active);
    })

    return strandLengthInputs;
  }

  private createSequenceBaseInput(): StringInput {
    const nucleotideBaseChoices = this.dataManager.fetchNucleotideBases();
    const defaultNucleotideBase = nucleotideBaseChoices[0];

    const sequenceBaseInput = ui.choiceInput('Sequence basis', defaultNucleotideBase, nucleotideBaseChoices, (value: string) => {
    });
    sequenceBaseInput.setTooltip('Nucleotide base to use for the sequence');
    return sequenceBaseInput;
  }

  private createPatternCommentInput(): StringInput {
    const patternCommentInput = ui.textInput('Comment', '', (value: string) => {});
    return patternCommentInput;
  }

  private createPatternSelectionBlock(): HTMLDivElement {
    const patternChoiceControls = new PatternChoiceControls(
      this.eventBus,
      this.dataManager,
    );
    return patternChoiceControls.getControlsContainer();
  }

  private createPatternNameInputBlock(): HTMLElement {
    const patternNameControls = new PatternNameControls(
      this.eventBus,
      this.patternConfiguration,
    );
    return patternNameControls.createPatternNameInputBlock();
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
    const userChoiceInput = this.createUserChoiceInput();
    const patternChoiceInput = this.getPatternChoiceInput();

    // todo: refactor this legacy solution
    patternChoiceInput.root.append(
      userChoiceInput.input,
      patternChoiceInput.input,
    );

    this.setPatternChoiceInputStyle(patternChoiceInput);

    const deletePatternButton = this.createDeletePatternButton();
    patternChoiceInput.addOptions(deletePatternButton);

    return patternChoiceInput;
  }

  private setPatternChoiceInputStyle(patternChoiceInput: StringInput): void {
    patternChoiceInput.setTooltip('Choose and apply pattern');
    patternChoiceInput.input.style.maxWidth = '120px';
    patternChoiceInput.input.style.marginLeft = '12px';
  }

  private createUserChoiceInput(): StringInput {
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

  private createDeletePatternButton(): HTMLDivElement {
    const button = ui.button(ui.iconFA('trash-alt'), () => this.eventBus.deletePattern(this.selectedPattern));

    return ui.div([ button ], 'ui-input-options');
  }

  private updatePatternChoiceInputContainer(): void {
    const patternInputs = this.getPatternInputs();
    $(this.patternChoiceContainer).empty();
    this.patternChoiceContainer.append(patternInputs.root);
  }
}

class PatternNameControls {
  constructor(
    private eventBus: EventBus,
    private patternConfiguration: PatternConfigurationManager,
  ) { }
  private patternName = 'Pattern';

  createPatternNameInputBlock(): HTMLElement {
    const patternNameInput = ui.textInput('Save as', this.patternName, (value: string) => this.handlePatternNameChange(value));

    this.handlePatternNameChange(patternNameInput.value);

    const savePatternButton = this.createSavePatternButton();

    patternNameInput.addOptions(savePatternButton);
    patternNameInput.setTooltip('Name of the pattern');
    return patternNameInput.root;
  }

  private createSavePatternButton(): HTMLElement {
    const savePatternButton = ui.bigButton('Save', () => this.processSaveButtonClick());
    return savePatternButton;
  }

  private handlePatternNameChange(patternName: string): void {
    this.patternName = patternName;
    this.patternConfiguration.setPatternName(this.patternName);
  }

  private processSaveButtonClick(): void {
    if (this.patternName === '') {
      grok.shell.warning(`Insert pattern name`);
      return;
    }
    grok.shell.info(`Pattern ${this.patternName} saved`);
    this.eventBus.requestPatternSave(this.patternName);
  }
}

class TableControlsManager {
  constructor(
    private eventBus: EventBus,
    private patternConfiguration: PatternConfigurationManager,
    private dataManager: PatternAppDataManager,
  ) {
    this.tableList = [];
    this.subscribeToTableEvents();
  }

  private tableList: DG.DataFrame[];
  private tableInputControlsContainer = ui.div([]);
  private columnControlsContainer = ui.div([]);

  private subscribeToTableEvents(): void {
    const observables = [
      grok.events.onTableAdded,
      grok.events.onTableRemoved,
    ];

    const handlers = [
      (table: DG.DataFrame) => this.addToTableList(table),
      (table: DG.DataFrame) => this.removeFromTableList(table),
    ];

    observables.forEach((observable, idx) => {
      const handle = handlers[idx];
      observable.subscribe((table: DG.DataFrame) => {
        console.log(`table ${table.name} added/removed`);
        handle(table);
        this.updateTableInputControls();
        // this.updateColumnControls(table);
      });
    });
  }

  private addToTableList(table: DG.DataFrame): void {
    this.tableList.push(table);
  }

  private removeFromTableList(table: DG.DataFrame): void {
    this.tableList = this.tableList.filter((item) => item.name !== table.name);
  }

  createUIComponents(): HTMLElement[] {
    const title = ui.h1('Convert');

    const tableInput = this.createTableInput();
    this.tableInputControlsContainer.append(tableInput);

    // this.updateColumnControls(this.tableList[0]);

    return [
      title,
      this.tableInputControlsContainer,
      this.columnControlsContainer,
    ];
  }

  private createTableInput(): HTMLElement {
    const tableInput = ui.tableInput('Tables', this.tableList[0], this.tableList, (table: DG.DataFrame) => {
      // this.updateColumnControls(table);
    });
    return tableInput.root;
  }

  private updateTableInputControls(): void {
    const newTableInput = this.createTableInput();
    $(this.tableInputControlsContainer).empty();
    this.tableInputControlsContainer.append(newTableInput);
  }

  // private createStrandColumnInput(table: DG.DataFrame): Record<StrandType, HTMLElement> {
  //   const columns = table ? table.columns.names() : [];
  //   console.log(`cols:`, columns);
  //   const strandColumnInput = Object.fromEntries(STRANDS.map((strand) => {
  //     const input = ui.choiceInput(`${STRAND_LABEL[strand]} column`, columns[0], columns, (colName: string) => { });
  //     return [strand, input.root];
  //   })) as Record<StrandType, HTMLElement>;
  //   return strandColumnInput;
  // }

  // private createIdColumnInput(table: DG.DataFrame): HTMLElement {
  //   const columns = table ? table.columns.names() : [];
  //   console.log(`cols:`, columns);
  //   const idColumnInput = ui.choiceInput('ID column', columns[0], columns, (colName: string) => { });
  //   return idColumnInput.root;
  // }

  // private updateColumnControls(table: DG.DataFrame): void {
  //   const strandColumnInput = this.createStrandColumnInput(table);
  //   const senseStrandColumnInput = strandColumnInput[SENSE_STRAND];
  //   const antisenseStrandColumnInput = strandColumnInput[ANTISENSE_STRAND];

  //   const idColumnInput = this.createIdColumnInput(table);

  //   $(this.columnControlsContainer).empty();
  //   this.columnControlsContainer.append(
  //     senseStrandColumnInput,
  //     antisenseStrandColumnInput,
  //     idColumnInput,
  //   );
  // }
}
