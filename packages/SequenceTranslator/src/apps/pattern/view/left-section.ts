/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';


import {MAX_SEQUENCE_LENGTH, OTHER_USERS, STRAND, STRANDS, STRAND_LABEL} from '../model/const';
import {StrandType, PatternNameExistsError, PatternExistsError} from '../model/types';

import {NumberInput, StringInput} from './types';

import $ from 'cash-dom';
import {PatternDefaultsProvider} from '../model/defaults-provider';
import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {PatternEditorDialog} from './pattern-editor';
import {isDialogOpen} from './pattern-editor';

export class PatternAppLeftSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
    private defaults: PatternDefaultsProvider
  ) { };

  getLayout(): HTMLDivElement {
    const patternControlsManager = new PatternControlsManager(
      this.eventBus,
      this.dataManager,
      this.defaults
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
    private defaultState: PatternDefaultsProvider
  ) { }

  createUIComponents(): HTMLElement[] {
    const antisenseStrandToggle = this.createAntisenseStrandToggle();
    const strandLengthInputs = this.createStrandLengthInputs();

    const senseStrandLengthInput = strandLengthInputs[STRAND.SENSE].root;
    const antisenseStrandLengthInput = strandLengthInputs[STRAND.ANTISENSE].root;

    const sequenceBaseInput = this.createSequenceBaseInput().root;
    const patternCommentInput = this.createPatternCommentInput().root;
    const patternSelectionBlock = this.createPatternSelectionBlock();
    const patternNameInputBlock = this.createPatternNameInputBlock();

    const editPatternButton = this.createEditPatternButton();
    // const savePatternButton = this.createSavePatternButton();

    return [
      ui.h1('Select'),
      patternSelectionBlock,

      ui.h1('Edit'),
      antisenseStrandToggle,
      senseStrandLengthInput,
      antisenseStrandLengthInput,
      sequenceBaseInput,
      patternNameInputBlock,
      patternCommentInput,
      ui.buttonsInput([
        editPatternButton,
        // savePatternButton,
      ]),
    ];
  }

  private createEditPatternButton(): HTMLButtonElement {
    const editPatternButton = ui.button(
      'Tweak',
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

    patternCommentInput.onInput(
      () => this.eventBus.updateComment(patternCommentInput.value!)
    );

    this.eventBus.patternLoaded$.subscribe(() => {
      patternCommentInput.value = this.eventBus.getComment();
    });

    return patternCommentInput;
  }

  private createPatternSelectionBlock(): HTMLDivElement {
    const patternChoiceControls = new PatternChoiceControls(
      this.eventBus,
      this.dataManager
    );
    return patternChoiceControls.getControlsContainer();
  }

  private createPatternNameInputBlock(): HTMLElement {
    const patternNameControls = new PatternNameControls(
      this.eventBus,
      this.dataManager
    );
    return patternNameControls.createPatternNameInputBlock();
  }
}

class PatternChoiceControls {
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

  getControlsContainer(): HTMLDivElement {
    const patternInputs = this.getPatternInputs();

    const container = ui.div([patternInputs.root]);
    this.eventBus.patternListUpdated$.subscribe(() => this.updatePatternChoiceInputContainer(container));

    return container;
  }

  private getPatternInputs(): StringInput {
    const userChoiceInput = this.createUserChoiceInput();
    const patternChoiceInput = this.createPatternChoiceInput();

    const patternInputContainer = ui.div([patternChoiceInput.input]);
    this.eventBus.userSelection$.subscribe(() => {
      $(patternInputContainer).empty();
      patternInputContainer.append(this.createPatternChoiceInput().input);
    });

    patternChoiceInput.root.append(
      userChoiceInput.input,
      patternInputContainer
    );

    this.setPatternChoiceInputStyle(patternChoiceInput);

    const deletePatternButton = this.createDeletePatternButton();
    this.eventBus.userSelection$.subscribe(() => {
      $(deletePatternButton).toggle(this.isCurrentUserSelected());
    });
    patternChoiceInput.addOptions(deletePatternButton);

    return patternChoiceInput;
  }

  private setPatternChoiceInputStyle(patternChoiceInput: StringInput): void {
    patternChoiceInput.setTooltip('Select and apply pattern');
    $(patternChoiceInput.input).css('max-width', '142px');
    patternChoiceInput.input.style.marginLeft = '12px';
  }

  private createUserChoiceInput(): StringInput {
    const currentUser = this.dataManager.getCurrentUserName();
    const possibleValues = [currentUser + ' (me)', OTHER_USERS];

    const userChoiceInput = ui.choiceInput(
      '',
      this.eventBus.getSelectedUser(),
      possibleValues,
      (userName: string) => this.eventBus.selectUser(userName)
    );
    this.setUserChoiceInputStyle(userChoiceInput);

    return userChoiceInput;
  }

  private setUserChoiceInputStyle(userChoiceInput: StringInput): void {
    userChoiceInput.setTooltip('Select pattern author');
    $(userChoiceInput.input).css('max-width', '140px');
  }

  private createPatternChoiceInput(): StringInput {
    const patternList = this.isCurrentUserSelected() ?
      [' '].concat(this.dataManager.getCurrentUserPatternNames()) :
      this.dataManager.getOtherUsersPatternNames();
    this.selectedPattern = patternList[0] || '<default>';
    this.eventBus.requestPatternLoad(this.selectedPattern);
    const choiceInput = ui.choiceInput(
      'Load pattern',
      this.selectedPattern,
      patternList
    );
    $(choiceInput.input).css({
      'max-width': '100px',
      'min-width': '100px',
    });
    choiceInput.onInput(
      () => this.eventBus.requestPatternLoad(choiceInput.value!)
    );
    return choiceInput;
  }

  private createDeletePatternButton(): HTMLDivElement {
    const button = ui.button(
      ui.iconFA('trash-alt'),
      () => this.showDeletePatternDialog()
    );

    ui.tooltip.bind(button, 'Delete pattern from user storage');

    return ui.div([button], 'ui-input-options');
  }

  private showDeletePatternDialog(): void {
    const dialog = ui.dialog('Delete pattern');
    dialog.add(ui.divText(`Are you sure you want to delete pattern ${this.selectedPattern}?`));
    dialog.onOK(() => this.eventBus.deletePattern(this.selectedPattern));
    dialog.show();
  }

  private updatePatternChoiceInputContainer(container: HTMLElement): void {
    const patternInputs = this.getPatternInputs();
    $(container).empty();
    container.append(patternInputs.root);
  }
}

class PatternNameControls {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager
  ) { }

  createPatternNameInputBlock(): HTMLElement {
    const patternNameInput = ui.textInput(
      'Pattern name',
      this.eventBus.getPatternName()
    );

    patternNameInput.onInput(
      () => this.eventBus.updatePatternName(patternNameInput.value)
    );
    this.eventBus.patternLoaded$.subscribe(() => {
      patternNameInput.value = this.eventBus.getPatternName();
    });

    // const savePatternButton = this.createSavePatternButton();

    // patternNameInput.addOptions(savePatternButton);
    // patternNameInput.setTooltip('Name of the pattern');
    return patternNameInput.root;
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
    // grok.shell.info(`Table ${table?.name} selected`);
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
      idColumnInput
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
    // console.log(`cols:`, columns);
    const idColumnInput = ui.choiceInput('ID column', columns[0], columns, (colName: string) => { });
    return idColumnInput.root;
  }
}
