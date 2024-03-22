/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';


import {MAX_SEQUENCE_LENGTH, OTHER_USERS, STRAND, STRANDS, STRAND_LABEL} from '../model/const';
import {StrandType} from '../model/types';

import './style.css';

import {NumberInput, StringInput} from './types';

import $ from 'cash-dom';
import {PatternDefaultsProvider} from '../model/defaults-provider';
import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {isDialogOpen, PatternEditorDialog} from './pattern-editor';

export class PatternAppLeftSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
    private defaults: PatternDefaultsProvider
  ) { };

  getLayout(): HTMLDivElement {
    const loadControlsManager = new PatternLoadControlsManager(
      this.eventBus,
      this.dataManager
    );

    const editControlsManager = new PatternEditControlsManager(
      this.eventBus,
      this.defaults
    );
    const tableControlsManager = new TableControlsManager(this.eventBus);

    const loadControls = loadControlsManager.createUIComponents();
    const editControls = editControlsManager.createUIComponents();
    const tableControls = tableControlsManager.createUIComponents();

    const loadControlsContainer = ui.div(loadControls);
    $(loadControlsContainer).css({'padding-bottom': '20px'});

    const form = ui.div([
      ...editControls,
      ...tableControls
    ], 'ui-form');

    const container = ui.div([
      loadControlsContainer,
      form
    ]);
    $(container).css({'padding': '25px'});

    const layout = ui.box(container, {style: {'maxWidth': '450px'}});
    return layout;
  }
}

class PatternEditControlsManager {
  constructor(
    private eventBus: EventBus,
    private defaultState: PatternDefaultsProvider
  ) { }

  createUIComponents(): HTMLElement[] {
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
          new PatternEditorDialog(this.eventBus).open();
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

    sequenceBaseInput.setTooltip('Most frequent nucleobase in the sequence');
    return sequenceBaseInput;
  }

  private createPatternCommentInput(): StringInput {
    const patternCommentInput = ui.textInput(
      'Comment',
      ''
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

class PatternLoadControlsManager {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager
  ) {
    this.eventBus.patternLoadRequested$.subscribe((value: string) => this.handlePatternChoice(value));

    const defaultUser = this.dataManager.getCurrentUserName();
    this.eventBus.selectUser(defaultUser);

    const defaultPattern = this.dataManager.getCurrentUserPatternNames()[0];
    this.selectedPattern = defaultPattern;
  }

  private selectedPattern: string;

  private async handlePatternChoice(patternName: string): Promise<void> {
    const patternConfiguration = await this.dataManager.getPatternConfig(patternName, this.isCurrentUserSelected());
    this.eventBus.setPatternConfig(patternConfiguration);
    this.selectedPattern = patternName;

    this.eventBus.updateControlsUponPatternLoaded();
  }

  private isCurrentUserSelected(): boolean {
    return this.eventBus.getSelectedUser() !== OTHER_USERS;
  }

  createUIComponents(): HTMLElement[] {
    const inputsContainer = this.getPatternInputsContainer();
    return [
      ui.h1('Load'),
      inputsContainer,
    ];
  }

  private getPatternInputsContainer(): HTMLDivElement {
    const inputsContainer = ui.divH(this.getPatternInputs());

    this.eventBus.patternListUpdated$.subscribe(() => {
      $(inputsContainer).empty();
      $(inputsContainer).append(this.getPatternInputs());
    });

    return inputsContainer;
  }

  private getPatternInputs(): HTMLElement[] {
    const userChoiceInput = this.createUserChoiceInput();
    const patternChoiceInputContainer = this.createPatternChoiceInputContainer();
    const deletePatternButton = this.createDeletePatternButton();

    return [
      userChoiceInput.root,
      patternChoiceInputContainer,
      deletePatternButton
    ];
  }

  private createPatternChoiceInputContainer(): HTMLDivElement {
    const patternChoiceInput = this.createPatternChoiceInput();
    const patternChoiceInputContainer = ui.div([patternChoiceInput.root]);
    this.eventBus.userSelection$.subscribe(() => {
      $(patternChoiceInputContainer).empty();
      $(patternChoiceInputContainer).append(this.createPatternChoiceInput().root);
    });

    return patternChoiceInputContainer;
  }

  private createUserChoiceInput(): StringInput {
    const currentUser = this.dataManager.getCurrentUserName();
    const possibleValues = [currentUser + ' (me)', OTHER_USERS];

    const userChoiceInput = ui.choiceInput(
      'Author',
      this.eventBus.getSelectedUser(),
      possibleValues,
      (userName: string) => this.eventBus.selectUser(userName)
    );
    this.setUserChoiceInputStyle(userChoiceInput);
    userChoiceInput.setTooltip('Select pattern author');

    return userChoiceInput;
  }

  private setUserChoiceInputStyle(userChoiceInput: StringInput): void {
    $(userChoiceInput.input).css({
      'max-width': '100px',
      'min-width': '100px',
    });
    $(userChoiceInput.root).css({
      'padding-right': '30px',
      'padding-left': '30px',
    });
  }

  private createPatternChoiceInput(): StringInput {
    const patternList = this.isCurrentUserSelected() ?
      [' '].concat(this.dataManager.getCurrentUserPatternNames()) :
      this.dataManager.getOtherUsersPatternNames();
    this.selectedPattern = patternList[0] || '<default>';
    this.eventBus.requestPatternLoad(this.selectedPattern);

    const choiceInput = ui.choiceInput('Pattern', this.selectedPattern, patternList);
    choiceInput.setTooltip('Select pattern to load');

    $(choiceInput.input).css({
      'max-width': '100px',
      'min-width': '100px',
    });

    choiceInput.onInput(
      () => this.eventBus.requestPatternLoad(choiceInput.value!)
    );
    return choiceInput;
  }

  private createDeletePatternButton(): HTMLButtonElement {
    const button = ui.button(
      ui.iconFA('trash-alt'),
      () => this.showDeletePatternDialog()
    );

    ui.tooltip.bind(button, 'Delete pattern from user storage');

    this.eventBus.userSelection$.subscribe(() => {
      $(button).toggle(this.isCurrentUserSelected());
    });

    return button;
  }

  private showDeletePatternDialog(): void {
    const dialog = ui.dialog('Delete pattern');
    dialog.add(ui.divText(`Are you sure you want to delete pattern ${this.selectedPattern}?`));
    dialog.onOK(() => this.eventBus.requestPatternDeletion(this.selectedPattern));
    dialog.show();
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
    const title = ui.h1('Bulk convert');
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
    if (this.availableTables.some((availableTable: DG.DataFrame) => availableTable.name === table.name))
      return;


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

    const idColumnInput = this.createIdColumnInput();

    $(this.columnControlsContainer).empty();
    this.columnControlsContainer.append(
      senseStrandColumnInput,
      antisenseStrandColumnInput,
      idColumnInput
    );
  }

  private createStrandColumnInput(): Record<StrandType, HTMLElement> {
    const columns = this.selectedTable ? this.selectedTable.columns.names() : [];
    const strandColumnInput = Object.fromEntries(STRANDS.map((strand) => {
      const input = ui.choiceInput(`${STRAND_LABEL[strand]} column`, columns[0], columns, () => { });
      return [strand, input.root];
    })) as Record<StrandType, HTMLElement>;
    return strandColumnInput;
  }

  private createIdColumnInput(): HTMLElement {
    const columns = this.selectedTable ? this.selectedTable.columns.names() : [];
    const idColumnInput = ui.choiceInput('ID column', columns[0], columns, () => { });
    return idColumnInput.root;
  }
}
