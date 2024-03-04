/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { STRAND, STRAND_LABEL, STRANDS, OTHER_USERS } from '../model/const';
import { StrandType } from '../model/types';

import {StringInput, NumberInput} from './types';

import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {PatternDefaultsProvider} from '../model/defaults-provider';
import $ from 'cash-dom';

export class PatternAppLeftSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
    private defaults: PatternDefaultsProvider,
  ) { };

  getLayout(): HTMLDivElement {
    const patternControlsManager = new PatternControlsManager(
      this.eventBus,
      this.dataManager,
      this.defaults,
    );
    const tableControlsManager = new TableControlsManager(this.eventBus);

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
    private dataManager: PatternAppDataManager,
    private defaultState: PatternDefaultsProvider,
  ) { }

  createUIComponents(): HTMLElement[] {
    const title = ui.h1('Pattern');

    const antisenseStrandToggle = this.createAntisenseStrandToggle();
    const strandLengthInputs = this.createStrandLengthInputs();

    const senseStrandLengthInput = strandLengthInputs[STRAND.SENSE].root;
    const antisenseStrandLengthInput = strandLengthInputs[STRAND.ANTISENSE].root;

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
      `${STRAND_LABEL[STRAND.ANTISENSE]} strand`,
      true,
      (isActive: boolean) => this.eventBus.toggleAntisenseStrand(isActive)
    );
    toggleAntisenseStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    return toggleAntisenseStrand.root;
  }

  private createStrandLengthInputs(): Record<string, NumberInput> {
    const createStrandLengthInput = (strand: StrandType) => {

      const sequenceLength = this.eventBus.getNucleotideSequences()[strand].length;

      const input = ui.intInput(
        `${STRAND_LABEL[strand]} length`,
        sequenceLength,
        (value: number) => updateStrandLengthInputs(strand, value, input)
      );
      input.setTooltip(`Length of ${STRAND_LABEL[strand].toLowerCase()}, including overhangs`);
      return [strand, input];
    }

    const updateStrandLengthInputs = (strand: StrandType, length: number, input: DG.InputBase) => {
      if (length < 0 || length > this.defaultState.getMaxSequenceLength()) {
        const oldValue = this.eventBus.getNucleotideSequences()[strand].length; 
        grok.shell.warning(`Sequence length must be between 0 and ${this.defaultState.getMaxSequenceLength()}`);
        input.value = oldValue;
        return;
      }

      this.eventBus.changeStrandLength(strand, length);
    }

    const strandLengthInputs = Object.fromEntries(
      STRANDS.map((strand) => createStrandLengthInput(strand))
    );

    this.eventBus.antisenseStrandToggled$.subscribe((active: boolean) => {
      $(strandLengthInputs[STRAND.ANTISENSE].root).toggle(active);
    })

    return strandLengthInputs;
  }

  private createSequenceBaseInput(): StringInput {
    const availableNucleoBases = this.defaultState.fetchAvailableNucleotideBases();
    const defaultNucleotideBase = this.defaultState.fetchDefaultNucleobase();

    const sequenceBaseInput = ui.choiceInput('Sequence basis', defaultNucleotideBase, availableNucleoBases, (base: string) => this.eventBus.changeSequenceBase(base));

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
      // this.patternConfiguration,
    );
    return patternNameControls.createPatternNameInputBlock();
  }
}

class PatternChoiceControls {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
  ) {
    this.eventBus.patternLoadRequested$.subscribe((value: string) => this.handlePatternChoice(value));

    const defaultUser = this.dataManager.getCurrentUserName();
    this.selectedUser = defaultUser;

    const defaultPattern = this.dataManager.getCurrentUserPatternNames()[0];
    this.selectedPattern = defaultPattern;

    this.patternChoiceContainer = ui.div([]);
    this.eventBus.patternListUpdated$.subscribe(() => this.updatePatternChoiceInputContainer());
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
    // private patternConfiguration: PatternConfigurationManager,
  ) { }

  createPatternNameInputBlock(): HTMLElement {
    const patternNameInput = ui.textInput(
      'Save as',
      this.eventBus.getPatternName(),
      (newValue: string) => this.eventBus.updatePatternName(newValue)
    );

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
    this.eventBus.updatePatternName(patternName);
  }

  private processSaveButtonClick(): void {
    if (this.eventBus.getPatternName() === '') {
      grok.shell.warning(`Insert pattern name`);
      return;
    }
    grok.shell.info(`Pattern ${this.eventBus.getPatternName()} saved`);
    this.eventBus.requestPatternSave();
  }
}

class TableControlsManager {
  private tableInputManager: TableInputManager;
  private columnInputManager: ColumnInputManager;

  constructor(eventBus: EventBus) {
    this.tableInputManager = new TableInputManager(eventBus);
    this.columnInputManager = new ColumnInputManager(eventBus);
  }

  createUIComponents(): HTMLElement[] {
    const title = ui.h1('Convert');
    const tableInput = this.tableInputManager.getTableInputContainer();
    const columnControls = this.columnInputManager.getColumnControlsContainer();

    const convertSequenceButton = ui.bigButton('Convert', () => this.processConvertButtonClick());

    return [
      title,
      tableInput,
      columnControls,
      ui.buttonsInput([
        convertSequenceButton,
      ]),
    ];
  }

  private processConvertButtonClick(): void {
  }
}

class TableInputManager {
  private availableTables: DG.DataFrame[] = [];
  private tableInputContainer: HTMLDivElement = ui.div([]);

  constructor(private eventBus: EventBus) {
    this.subscribeToTableEvents();
    this.refreshTableInput();
  }

  getTableInputContainer(): HTMLDivElement {
    return this.tableInputContainer;
  }

  private subscribeToTableEvents(): void {
    grok.events.onTableAdded.subscribe((table: DG.DataFrame) => this.handleTableAdded(table));
    grok.events.onTableRemoved.subscribe((table: DG.DataFrame) => this.handleTableRemoved(table));
    this.eventBus.tableSelectionChanged$.subscribe(() => this.handleTableChoice());
  }

  private handleTableAdded(table: DG.DataFrame): void {
    if (this.availableTables.some((availableTable: DG.DataFrame) => availableTable.name === table.name)) {
      return;
    }

    this.availableTables.push(table);

    this.refreshTableInput();
  }

  private handleTableRemoved(removedTable: DG.DataFrame): void {
    this.availableTables = this.availableTables.filter((table: DG.DataFrame) => table.name !== removedTable.name);

    this.refreshTableInput();
  }

  private refreshTableInput(): void {
    const tableInput = this.createTableInput();
    $(this.tableInputContainer).empty();
    this.tableInputContainer.append(tableInput.root);
  }

  private createTableInput(): DG.InputBase<DG.DataFrame | null> {
    const currentSelection = this.eventBus.getTableSelection();

    const tableInput = ui.tableInput(
      'Tables',
      currentSelection,
      this.availableTables,
      (table: DG.DataFrame) => this.eventBus.selectTable(table));
    return tableInput;
  }

  private handleTableChoice(): void {
    const table = this.eventBus.getTableSelection();
    if (!table) return;
    if (!this.isTableDisplayed(table))
      this.displayTable(table);
    grok.shell.info(`Table ${table?.name} selected`);
  }

  private isTableDisplayed(table: DG.DataFrame): boolean {
    return grok.shell.tableNames.includes(table.name);
  }

  private displayTable(table: DG.DataFrame): void {
    const previousView = grok.shell.v;
    grok.shell.addTableView(table);
    grok.shell.v = previousView;
  }
}

class ColumnInputManager {
  private columnControlsContainer: HTMLDivElement = ui.div([]);

  constructor(private eventBus: EventBus) {
    this.eventBus.tableSelectionChanged$.subscribe(() => this.handleTableChoice());
    this.refreshColumnControls();
  }

  getColumnControlsContainer(): HTMLDivElement {
    return this.columnControlsContainer;
  }

  private get selectedTable(): DG.DataFrame | null {
    return this.eventBus.getTableSelection();
  }

  private handleTableChoice(): void {
    this.refreshColumnControls();
    grok.shell.info(`Table ${this.selectedTable?.name} selection from column input manager`);
  }

  private refreshColumnControls(): void {
    const strandColumnInput = this.createStrandColumnInput();
    const senseStrandColumnInput = strandColumnInput[STRAND.SENSE];
    const antisenseStrandColumnInput = strandColumnInput[STRAND.ANTISENSE];
    // $(antisenseStrandColumnInput).toggle();

    const idColumnInput = this.createIdColumnInput();

    $(this.columnControlsContainer).empty();
    this.columnControlsContainer.append(
      senseStrandColumnInput,
      antisenseStrandColumnInput,
      idColumnInput,
    );
  }

  private createStrandColumnInput(): Record<StrandType, HTMLElement> {
    const columns = this.selectedTable ? this.selectedTable.columns.names() : [];
    const strandColumnInput = Object.fromEntries(STRANDS.map((strand) => {
      const input = ui.choiceInput(`${STRAND_LABEL[strand]} column`, columns[0], columns, (colName: string) => { });
      return [strand, input.root];
    })) as Record<StrandType, HTMLElement>;
    return strandColumnInput;
  }

  private createIdColumnInput(): HTMLElement {
    const columns = this.selectedTable ? this.selectedTable.columns.names() : [];
    console.log(`cols:`, columns);
    const idColumnInput = ui.choiceInput('ID column', columns[0], columns, (colName: string) => { });
    return idColumnInput.root;
  }
}
