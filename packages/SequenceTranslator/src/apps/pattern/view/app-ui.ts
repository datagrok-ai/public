/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {AXOLABS_STYLE_MAP} from '../../common/data-loader/json-loader';
import {
  DEFAULT_PATTERN_CONFIG as DEFAULT, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, STRAND, STRAND_LABEL, STRANDS, TERMINI, PATTERN_KEY, TERMINUS
} from '../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../model/types';
import {isOverhangNucleotide} from '../model/helpers';
import {generateExample, translateSequence, getShortName, isPatternCreatedByCurrentUser, findDuplicates, addColumnWithIds, addColumnWithTranslatedSequences} from '../model/oligo-pattern';
// import {renderNucleotidePattern} from './svg-utils/render-svg';

import { BooleanInput, StringInput} from './types';
import {LegacyPatternConfig} from '../model/types';
// todo: remove ts-ignore
//@ts-ignore
import * as svgExport from 'save-svg-as-png';
import $ from 'cash-dom';

type Bases = Record<typeof STRANDS[number], string>;
type PtoLinkageActive = Record<typeof STRANDS[number], boolean[]>;
type TerminalModifications = Record<typeof STRANDS[number], Record<typeof TERMINI[number], string>>;

const enum MODIFICATION_CATEGORY {
  PTO,
  BASIS
}

export class PatternLayoutHandler {
  get htmlDivElement() {
    function updateStrandModificationControls(strand: string) {
      clearModificationControls(strand);
      extendStrandInputsForNewSequenceLength(strand);
      let nucleotideCounter = 0;

      for (let i = 0; i < getStrandLength(strand); i++) {
        updateSinglePhosphorothioateLinkageInput(strand, i);
        updateSingleNucleotideBaseInput(strand, i);

        if (!isOverhangNucleotide(getBaseInputValue(strand, i))) {
          nucleotideCounter++;
        }

        const modificationControlGroup = generateSingleModificationControlGroup(strand, nucleotideCounter, i);
        modificationControls[strand].append(modificationControlGroup);
      }
    }

    function clearModificationControls(strand: string) {
      modificationControls[strand].innerHTML = '';
    }

    function extendStrandInputsForNewSequenceLength(strand: string) {
      const extensionLength = maxStrandLength[strand] - nucleobaseInputs[strand].length;
      ptoLinkageInputs[strand] = ptoLinkageInputs[strand].concat(Array(extensionLength).fill(isFullyPtoInput));
      nucleobaseInputs[strand] = nucleobaseInputs[strand].concat(Array(extensionLength).fill(sequenceBase));
    }

    function getStrandLength(strand: string) {
      return strandLengthInputs[strand].value!;
    }

    function updateSinglePhosphorothioateLinkageInput(strand: string, index: number) {
      ptoLinkageInputs[strand][index] = ui.boolInput('', ptoLinkageInputs[strand][index].value!, () => {
        refreshSvgDisplay();
        refreshOutputExamples();
      });
    }

    function updateSingleNucleotideBaseInput(strand: string, index: number) {
      nucleobaseInputs[strand][index] = ui.choiceInput('', getBaseInputValue(strand, index), nucleotideBaseChoices, (value: string) => {
        handleBaseInputChange(value);
        updateStrandModificationControls(STRAND.ANTISENSE);
        refreshSvgDisplay();
        refreshOutputExamples();
      });

      $(nucleobaseInputs[strand][index].root).addClass('st-pattern-choice-input');
    }

    function getBaseInputValue(strand: string, index: number) {
      return nucleobaseInputs[strand][index].value!;
    }

    function handleBaseInputChange(value: string) {
      if (!modificationsWithNumericLabels.includes(value)) {
        toggleModificationInList(value);
      }
    }

    function updateNumberedModificationsList(value: string, shouldAdd: boolean): void {
      const index = modificationsWithNumericLabels.indexOf(value);
      if (shouldAdd && index === -1) {
        modificationsWithNumericLabels.push(value);
      } else if (!shouldAdd && index > -1) {
        modificationsWithNumericLabels.splice(index, 1);
      }
    }

    function appendToNumberedModificationsContainer(
      value: string,
      onChangeCallback: (value: string, shouldAdd: boolean) => void
    ): void {
      const boolInput = ui.boolInput(value, true, (boolV: boolean) => onChangeCallback(value, boolV));
      numericLabelTogglesContainer.append(
        ui.divText('', {style: {width: '25px'}}),
        boolInput.root,
      );
    }

    function toggleModificationInList(value: string): void {
      const shouldAdd = !modificationsWithNumericLabels.includes(value);
      updateNumberedModificationsList(value, shouldAdd);
      if (shouldAdd) {
        appendToNumberedModificationsContainer(value, (val, shouldAdd) => {
          updateNumberedModificationsList(val, shouldAdd);
          refreshSvgDisplay();
        });
      }
    }

    function generateSingleModificationControlGroup(strand: string, nucleotideCounter: number, index: number) {
      const labelText = isOverhangNucleotide(getBaseInputValue(strand, index)) ? '' : String(nucleotideCounter);
      const labelUI = ui.div([ui.label(labelText)], {style: {width: '20px'}});
      const baseInputUI = ui.block75([nucleobaseInputs[strand][index].root]);
      const ptoLinkageUI = ui.div([ptoLinkageInputs[strand][index]]);

      return ui.divH([labelUI, baseInputUI, ptoLinkageUI], {style: {alignItems: 'center'}});
    }

    function refreshUIForNewSequenceLength() {
      if (checkSequenceLengthValidity()) {
        extendStrandsToNewLength();
        refreshUIComponents();
      } else {
        displayOutOfRangeDialog();
      }
    }

    function checkSequenceLengthValidity() {
      return Object.values(strandLengthInputs).every(input => input.value! < MAX_SEQUENCE_LENGTH);
    }

    function extendStrandsToNewLength() {
      STRANDS.forEach(strand => {
        if (strandLengthInputs[strand].value! > maxStrandLength[strand]) {
          maxStrandLength[strand] = strandLengthInputs[strand].value!;
        }
        updateStrandModificationControls(strand);
      });
    }

    function refreshUIComponents() {
      refreshSvgDisplay();
      refreshInputExamples();
      refreshOutputExamples();
    }

    function displayOutOfRangeDialog() {
      ui.dialog('Out of range')
        .add(ui.divText(`Sequence length should be less than ${MAX_SEQUENCE_LENGTH} due to UI constraints.`))
        .onOK(() => resetAllStrandLengthsToMax())
        .onCancel(() => resetAllStrandLengthsToMax())
        .showModal(false);
    }

    function resetAllStrandLengthsToMax() {
      Object.values(strandLengthInputs).forEach(input => input.value = MAX_SEQUENCE_LENGTH);
    }

    function updateInputValues(category: MODIFICATION_CATEGORY, newValue: boolean | string): void {
      const inputs = category === MODIFICATION_CATEGORY.PTO ? ptoLinkageInputs : nucleobaseInputs;

      for (let i = 0; i < STRANDS.length; i++) {
        const strand = STRANDS[i];
        // WARNING: replacing this with for (const ...) or .forEach() leads to a bug: some values are not updated !!!
        for (let j = 0; j < inputs[strand].length; j++) {
          const item = inputs[strand][j];
          item.value = newValue;
        }
      }
      refreshSvgDisplay();
    }

    function refreshInputExamples(): void {
      STRANDS.forEach((strand) => {
        if (strandColumnInput[strand].value === '') {
          inputExample[strand].value = generateExample(strandLengthInputs[strand].value!, sequenceBase.value!);
        }
      });
    }

    function refreshOutputExamples(): void {
      const conditions = [true, createAsStrand.value];
      STRANDS.forEach((strand, idx) => {
        if (conditions[idx]) {
          outputExample[strand].value = translateSequence(
            inputExample[strand].value,
            nucleobaseInputs[strand],
            ptoLinkageInputs[strand],
            terminalModificationInputs[strand][TERMINUS.FIVE_PRIME],
            terminalModificationInputs[strand][TERMINUS.THREE_PRIME],
            firstPto[strand].value!
          );
        }
      });
    }

    function getBaseInputValues(strand: string) {
      return nucleobaseInputs[strand].slice(0, strandLengthInputs[strand].value!).map(e => e.value!);
    }

    function getPtoLinkageValues(strand: string): boolean[] {
      return [firstPto[strand].value!, ...ptoLinkageInputs[strand].slice(0, strandLengthInputs[strand].value!).map(e => e.value!)];
    }

    function refreshSvgDisplay() {
      svgDisplayDiv.innerHTML = '';
      const baseInputValues = Object.fromEntries(STRANDS.map((strand) => [strand, getBaseInputValues(strand)]));
      const ptoLinkageValues = Object.fromEntries(STRANDS.map((strand) => [strand, getPtoLinkageValues(strand)]));

      const terminalModificationValues = Object.fromEntries(STRANDS.map((strand) => [strand, {
        [TERMINUS.THREE_PRIME]: terminalModificationInputs[strand][TERMINUS.THREE_PRIME].value!,
        [TERMINUS.FIVE_PRIME]: terminalModificationInputs[strand][TERMINUS.FIVE_PRIME].value!,
      }]));

      const patternSettings = {
        patternName: getShortName(patternNameInput.value),
        isAntisenseStrandIncluded: createAsStrand.value!,
        nucleotideSequences: baseInputValues,
        phosphorothioateLinkageFlags: ptoLinkageValues,
        strandTerminusModifications: terminalModificationValues,
        patternComment: patternCommentInput.value,
        nucleotidesWithNumericLabels: modificationsWithNumericLabels,
      } as PatternConfiguration;

      // const svg = renderNucleotidePattern(patternSettings);

      // svgDisplayDiv.append(ui.span([svg]));
    }

    function detectMostFrequentBase(array: string[]): string {
      const countMap: {[element: string]: number} = {};
      let mostFrequentElement = array[0];
      let highestCount = 1;

      array.forEach(element => {
        countMap[element] = (countMap[element] || 0) + 1;
        if (countMap[element] > highestCount) {
          mostFrequentElement = element;
          highestCount = countMap[element];
        }
      });

      return mostFrequentElement;
    }

    async function fetchPatternFromStorage(newName: string): Promise<LegacyPatternConfig> {
      const entities = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);
      return JSON.parse(entities[newName]);
    }

    async function applyPatternDataToUI(newName: string, patternConfig: LegacyPatternConfig): Promise<void> {
      sequenceBase.value = detectMostFrequentBase([...patternConfig[PATTERN_KEY.AS_BASES], ...patternConfig[PATTERN_KEY.SS_BASES]]);
      createAsStrand.value = (patternConfig[PATTERN_KEY.AS_BASES].length > 0);
      patternNameInput.value = newName;

      refreshBaseInputsFromPattern(patternConfig);
      refreshPtoLinkagesFromPattern(patternConfig);
      adjustStrandLengthsFromPattern(patternConfig);
      refreshTerminalModificationsFromPattern(patternConfig);

      patternCommentInput.value = patternConfig[PATTERN_KEY.COMMENT];
    }

    async function fetchAndUpdatePatternInUI(newName: string): Promise<void> {
      const progressIndicator = DG.TaskBarProgressIndicator.create('Loading pattern...');
      try {
        const patternObj = await fetchPatternFromStorage(newName);
        await applyPatternDataToUI(newName, patternObj);
      } catch (error) {
        console.error("Error parsing pattern and updating UI: ", error);
      } finally {
        progressIndicator.close();
      }
    }

    function refreshBaseInputsFromPattern(obj: LegacyPatternConfig): void {
      const fields = [PATTERN_KEY.SS_BASES, PATTERN_KEY.AS_BASES];
      STRANDS.forEach((strand, i) => {
        nucleobaseInputs[strand] = (obj[fields[i]] as string[]).map((base: string) => ui.choiceInput('', base, nucleotideBaseChoices));
      });
    }

    function refreshPtoLinkagesFromPattern(obj: LegacyPatternConfig): void {
      const fields = [PATTERN_KEY.SS_PTO, PATTERN_KEY.AS_PTO];
      STRANDS.forEach((strand, i) => {
        const ptoValues = obj[fields[i]] as boolean[];
        firstPto[strand].value = ptoValues[0];
        ptoLinkageInputs[strand] = ptoValues.slice(1).map((value: boolean) => ui.boolInput('', value));
      });
    }

    function adjustStrandLengthsFromPattern(obj: LegacyPatternConfig): void {
      const fields = [PATTERN_KEY.SS_BASES, PATTERN_KEY.AS_BASES];
      STRANDS.forEach((strand, i) => {
        strandLengthInputs[strand].value = obj[fields[i]].length;
      });
    }

    function refreshTerminalModificationsFromPattern(obj: LegacyPatternConfig): void {
      const field = [[PATTERN_KEY.SS_3, PATTERN_KEY.SS_5], [PATTERN_KEY.AS_3, PATTERN_KEY.AS_5]];
      STRANDS.forEach((strand, i) => {
        TERMINI.forEach((terminal, j) => {
          terminalModificationInputs[strand][terminal].value = obj[field[i][j]] as string;
        });
      });
    }

    function verifyUniformColumnLengths(colName: string): boolean {
      const col = tableInput.value!.getCol(colName);
      const areLengthsUniform = col.toList().every((value, index, array) =>
        index === 0 || value.length === array[index - 1].length || value.length === 0);

      if (!areLengthsUniform) {
        displayLengthMismatchDialog(colName);
      } else if (col.get(0).length !== strandLengthInputs[STRAND.SENSE].value) {
        displayLengthUpdatedDialog();
      }

      return areLengthsUniform;
    }

    function displayLengthMismatchDialog(colName: string) {
      const dialog = ui.dialog('Sequences lengths mismatch');
      $(dialog.getButton('OK')).hide();

      dialog
        .add(ui.divText('The sequence length should match the number of Raw sequences in the input file'))
        .add(ui.divText('\'ADD COLUMN\' to see sequences lengths'))
        .addButton('ADD COLUMN', () => addColumnWithSequenceLengths(colName, dialog))
        .show();
    }

    function addColumnWithSequenceLengths(colName: string, dialog: any) {
      const table = tableInput.value!;
      table.columns.addNewInt('Sequences lengths in ' + colName).init((j: number) => table.getCol(colName).get(j).length);
      grok.shell.info('Column with lengths added to \'' + table.name + '\'');
      dialog.close();
      grok.shell.v = grok.shell.getTableView(table.name);
    }

    function displayLengthUpdatedDialog() {
      ui.dialog('Length was updated by value from imported file')
        .add(ui.divText('Latest modifications may not take effect during translation'))
        .onOK(() => grok.shell.info('Lengths changed'))
        .show();
    }

    async function fetchCurrentUserName(): Promise<string> {
      const user = await grok.dapi.users.current();
      return ` (created by ${user.friendlyName})`;
    }

    async function savePatternAndNotifyUser(): Promise<void> {
      const currUserName = await fetchCurrentUserName();
      patternNameInput.value = createPatternNameForSaving(currUserName);

      const patternData = assemblePatternData();
      await grok.dapi.userDataStorage.postValue(USER_STORAGE_KEY, patternNameInput.value, JSON.stringify(patternData), false);

      grok.shell.info(`Pattern '${patternNameInput.value}' was successfully uploaded!`);
    }

    function createPatternNameForSaving(currUserName: string): string {
      return patternNameInput.stringValue.includes('(created by ') ?
        `${getShortName(patternNameInput.value)}${currUserName}` :
        `${patternNameInput.stringValue}${currUserName}`;
    }

    function assemblePatternData(): LegacyPatternConfig {
      const createBasesArray = (strand: string) => nucleobaseInputs[strand].slice(0, strandLengthInputs[strand].value!).map(e => e.value) as string[];
      const createPtoArray = (strand: string) => [firstPto[strand].value, ...ptoLinkageInputs[strand].slice(0, strandLengthInputs[strand].value!).map(e => e.value)] as boolean[];

      return {
        [PATTERN_KEY.SS_BASES]: createBasesArray(STRAND.SENSE),
        [PATTERN_KEY.AS_BASES]: createBasesArray(STRAND.ANTISENSE),
        [PATTERN_KEY.SS_PTO]: createPtoArray(STRAND.SENSE),
        [PATTERN_KEY.AS_PTO]: createPtoArray(STRAND.ANTISENSE),
        [PATTERN_KEY.SS_3]: terminalModificationInputs[STRAND.SENSE][TERMINUS.THREE_PRIME].value!,
        [PATTERN_KEY.SS_5]: terminalModificationInputs[STRAND.SENSE][TERMINUS.FIVE_PRIME].value!,
        [PATTERN_KEY.AS_3]: terminalModificationInputs[STRAND.ANTISENSE][TERMINUS.THREE_PRIME].value!,
        [PATTERN_KEY.AS_5]: terminalModificationInputs[STRAND.ANTISENSE][TERMINUS.FIVE_PRIME].value!,
        [PATTERN_KEY.COMMENT]: patternCommentInput.value,
      };
    }

    async function sortPatternsByUserOwnership(patternsData: LegacyPatternConfig): Promise<{ ownPatterns: string[], otherUsersPatterns: string[] }> {
      const ownPatterns: string[] = [];
      const otherUsersPatterns: string[] = [];

      for (const patternName of Object.keys(patternsData)) {
        if (await isPatternCreatedByCurrentUser(patternName))
          otherUsersPatterns.push(patternName);
        else
          ownPatterns.push(patternName);
      }

      return { ownPatterns, otherUsersPatterns };
    }

    function initializeLoadPatternInterface(loadPattern: StringInput, loadPatternDiv: HTMLElement, patternListChoiceInput: StringInput) {
      clearAndRepopulateUIElements(loadPatternDiv, loadPattern, patternListChoiceInput);
      styleLoadPatternInputField(loadPattern.input);
      loadPattern.setTooltip('Apply Existing Pattern');
      const deleteButton = createPatternDeletionButton(loadPattern);
      loadPattern.root.append(deleteButton);
    }

    function clearAndRepopulateUIElements(loadPatternDiv: HTMLElement, loadPattern: StringInput, patternListChoiceInput: StringInput) {
      loadPatternDiv.innerHTML = '';
      loadPattern.root.append(patternListChoiceInput.input, loadPattern.input);
      loadPatternDiv.append(loadPattern.root);
    }

    function styleLoadPatternInputField(loadPatternInput: HTMLElement) {
      loadPatternInput.style.maxWidth = '120px';
      loadPatternInput.style.marginLeft = '12px';
    }

    function createPatternDeletionButton(loadPattern: StringInput): HTMLElement {
      return ui.div([
        ui.button(ui.iconFA('trash-alt'), async () => {
          if (!loadPattern.value) {
            grok.shell.warning('Choose pattern to delete');
          } else if (await isPatternCreatedByCurrentUser(patternNameInput.value)) {
            grok.shell.warning('Cannot delete pattern, created by other user');
          } else {
            await removePatternFromStorage(loadPattern.value);
            await refreshPatternsList();
          }
        })
      ], 'ui-input-options');
    }

    async function removePatternFromStorage(patternName: string) {
      await grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, patternName, false);
      grok.shell.info(`Pattern '${patternName}' deleted`);
    }

    async function refreshPatternsList() {
      const patternsData = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);
      const { ownPatterns, otherUsersPatterns } = await sortPatternsByUserOwnership(patternsData);

      const currentUserName = (await grok.dapi.users.current()).friendlyName;
      const otherUsers = 'Other users';

      let patternChoiceInput = ui.choiceInput('Load pattern', '', ownPatterns, (value: string) => fetchAndUpdatePatternInUI(value));

      const userChoiceInput = ui.choiceInput(
        '', currentUserName, [currentUserName, otherUsers],
        (value: string) => patternListChoiceInputOnChange(value)
      );

      userChoiceInput.input.style.maxWidth = '142px';
      initializeLoadPatternInterface(patternChoiceInput, loadPatternDiv, userChoiceInput);

      function patternListChoiceInputOnChange(value: string) {
        const currentList = value === currentUserName ? ownPatterns : otherUsersPatterns;
        patternChoiceInput = ui.choiceInput('Load pattern', '', currentList, (value: string) => fetchAndUpdatePatternInUI(value));
        initializeLoadPatternInterface(patternChoiceInput, loadPatternDiv, userChoiceInput);
      }
    }

    async function checkIfPatternExistsInStorage(patternData: LegacyPatternConfig, patternName: string) {
      return Object.keys(patternData).includes(patternName);
    }

    async function displayPatternReplaceConfirmationDialog(patternName: string) {
      const dialog = ui.dialog('Pattern already exists');
      $(dialog.getButton('OK')).hide();
      dialog
        .add(ui.divText(`Pattern name '${patternName}' already exists.`))
        .add(ui.divText('Replace pattern?'))
        .addButton('YES', async () => {
          await deletePatternFromStorage(patternName);
          await savePatternAndNotifyUser();
          dialog.close();
        })
        .show();
    }

    async function deletePatternFromStorage(patternName: string) {
      return grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, patternName, false);
    }

    async function savePatternInUserDataStorage(): Promise<void> {
      const patternData = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);

      const patternName = patternNameInput.value;

      if (await checkIfPatternExistsInStorage(patternData, patternName)) {
        await displayPatternReplaceConfirmationDialog(patternName);
      } else {
        await savePatternAndNotifyUser();
      }

      await refreshPatternsList();
    }

    function validateStrandColumn(colName: string, strand: string): void {
      if (!verifyUniformColumnLengths(colName)) {
        return;
      }

      const firstSequence = fetchFirstSequenceFromColumn(colName);
      adjustStrandLengthIfRequired(firstSequence.length, strand);
      refreshInputExampleValue(firstSequence, strand);
    }

    function fetchFirstSequenceFromColumn(colName: string): string {
      return tableInput.value!.getCol(colName).get(0);
    }

    function adjustStrandLengthIfRequired(sequenceLength: number, strand: string): void {
      const currentStrandLength = strandLengthInputs[strand].value;
      if (sequenceLength !== currentStrandLength) {
        strandLengthInputs[strand].value = sequenceLength;
      }
    }

    function refreshInputExampleValue(sequence: string, strand: string): void {
      inputExample[strand].value = sequence;
    }

    function validateIdsColumn(colName: string) {
      const col = retrieveColumnFromTable(colName);
      checkColumnTypeForIntType(col);

      if (checkForDuplicateValues(col)) {
        showDuplicatesDialog(col);
      }
    }

    function retrieveColumnFromTable(colName: string): DG.Column {
      return tableInput.value!.getCol(colName);
    }

    function checkColumnTypeForIntType(col: DG.Column): void {
      if (col.type !== DG.TYPE.INT) {
        grok.shell.error('Column should contain integers only');
        // throw new Error('Invalid column type');
      }
    }

    function checkForDuplicateValues(col: DG.Column): boolean {
      const uniqueValues = col.categories.filter((e) => e !== '').length;
      const totalValues = col.toList().filter((e) => e !== '').length;
      return uniqueValues < totalValues;
    }

    function showDuplicatesDialog(col: DG.Column): void {
      const duplicates = findDuplicates(col.getRawData());
      const tableName = tableInput.value!.name;

      ui.dialog('Non-unique IDs')
      .add(ui.divText('Press \'OK\' to select rows with non-unique values'))
      .onOK(() => selectDuplicateRows(duplicates, col))
      .show();

      grok.shell.info(`Rows are selected in table '${tableName}'`);
    }

    function selectDuplicateRows(duplicates: number[], col: DG.Column): void {
      const selection = tableInput.value!.selection;
      selection.init((i: number) => duplicates.includes(col.get(i)));
      grok.shell.v = grok.shell.getTableView(tableInput.value!.name);
    }


    const nucleotideBaseChoices: string[] = Object.keys(AXOLABS_STYLE_MAP);
    const defaultNucleotideBase: string = nucleotideBaseChoices[0];
    const modificationsWithNumericLabels = [defaultNucleotideBase];
    const sequenceBase = ui.choiceInput('Sequence basis', defaultNucleotideBase, nucleotideBaseChoices, (value: string) => {
      // updateBases(value);
      updateInputValues(MODIFICATION_CATEGORY.BASIS, value);
      refreshOutputExamples();
    });

    function createFullyPtoInput(): BooleanInput {
      const fullyPto = ui.boolInput('Fully PTO', DEFAULT.PHOSPHOROTHIOATE, handleFullPtoInputChange);
      styleFullyPtoInputLabel(fullyPto.captionLabel);
      return fullyPto;
    }

    function styleFullyPtoInputLabel(label: HTMLElement): void {
      label.classList.add('ui-label-right');
      Object.assign(label.style, {
        textAlign: 'left',
        maxWidth: '100px',
        minWidth: '40px',
        width: 'auto'
      });
    }

    function handleFullPtoInputChange(value: boolean): void {
      STRANDS.forEach((strand) => { firstPto[strand].value = value; });
      // updatePto(value);
      updateInputValues(MODIFICATION_CATEGORY.PTO, value);
      refreshOutputExamples();
    }

    const isFullyPtoInput = createFullyPtoInput();

    const maxStrandLength = Object.fromEntries(STRANDS.map(
      (strand) => [strand, DEFAULT.SEQUENCE_LENGTH]
    ));
    const modificationControls = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.div([])]
    ));
    const ptoLinkageInputs = Object.fromEntries(STRANDS.map(
      (strand) => [strand, Array<BooleanInput>(DEFAULT.SEQUENCE_LENGTH)
        .fill(ui.boolInput('', DEFAULT.PHOSPHOROTHIOATE))]
    ));
    const nucleobaseInputs = Object.fromEntries(STRANDS.map(
      (strand) => {
        const choiceInputs = Array<StringInput>(DEFAULT.SEQUENCE_LENGTH)
          .fill(ui.choiceInput('', defaultNucleotideBase, nucleotideBaseChoices));
        return [strand, choiceInputs];
      }
    ));
    const strandLengthInputs = Object.fromEntries(STRANDS.map(
      (strand) => {
        const input = ui.intInput(`${STRAND_LABEL[strand]} length`, DEFAULT.SEQUENCE_LENGTH, () => refreshUIForNewSequenceLength());
        input.setTooltip(`Length of ${STRAND_LABEL[strand].toLowerCase()}, including overhangs`);
        return [strand, input];
      }));
    const selectedStrandColumn = Object.fromEntries(STRANDS.map((strand) => [strand, '']));
    const inputExample = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.textInput(
        ``, generateExample(strandLengthInputs[strand].value!, sequenceBase.value!))
      ]));

    const strandColumnInput = Object.fromEntries(STRANDS.map((strand) => {
      const input = ui.choiceInput(`${STRAND_LABEL[strand]} column`, '', [], (colName: string) => {
        validateStrandColumn(colName, strand);
        selectedStrandColumn[strand] = colName;
      });
      return [strand, input];
    }));

    function generateInitialPtoInputs(fullyPto: BooleanInput): Record<string,BooleanInput> {
      return Object.fromEntries(STRANDS.map((strand) => [strand, generateStrandSpecificPtoInput(strand, fullyPto)]));
    }

    function generateStrandSpecificPtoInput(strand: StrandType, fullyPto: BooleanInput): BooleanInput {
      const input = ui.boolInput(`First ${strand} PTO`, fullyPto.value!, refreshSvgDisplay);
      input.setTooltip(`ps linkage before first nucleotide of ${STRAND_LABEL[strand].toLowerCase()}`);
      styleStrandPtoInputLabel(input.captionLabel);
      return input;
    }

    function styleStrandPtoInputLabel(label: HTMLElement): void {
      label.classList.add('ui-label-right');
      Object.assign(label.style, {
        textAlign: 'left',
        maxWidth: '100px',
        minWidth: '40px',
        width: 'auto'
      });
    }

    const firstPto = generateInitialPtoInputs(isFullyPtoInput);

    function createTerminalModificationInputs() {
      return Object.fromEntries(STRANDS.map(strand => [strand, createStrandInputs(strand)]));
    }

    function createStrandInputs(strand: StrandType) {
      return Object.fromEntries(TERMINI.map((terminus: TerminalType) => [terminus, createTerminalModificationInput(strand, terminus)]));
    }

    function createTerminalModificationInput(strand: StrandType, terminus: TerminalType): StringInput {
      const label = `${strand} ${terminus} Modification`;
      const tooltip = `Additional ${strand} ${terminus} Modification`;
      const input = ui.stringInput(label, '', handleTerminalModificationInputChange);

      input.setTooltip(tooltip);
      return input;
    }

    function handleTerminalModificationInputChange(): void {
      refreshSvgDisplay();
      refreshOutputExamples();
    }

    const terminalModificationInputs = createTerminalModificationInputs();

    function createOutputExample(strand: StrandType) {
      const translatedSequence = translateSequence(
        inputExample[strand].value,
        nucleobaseInputs[strand],
        ptoLinkageInputs[strand],
        terminalModificationInputs[strand][TERMINUS.THREE_PRIME],
        terminalModificationInputs[strand][TERMINUS.FIVE_PRIME],
        firstPto[strand].value!
      );

      const input = ui.textInput('', translatedSequence);

      Object.assign(input.input.style, {
        minWidth: 'none',
        flexGrow: '1',
      });
      $(input.root.lastChild).css('height', 'auto');

      return input;
    }

    const outputExample = Object.fromEntries(
      STRANDS.map((strand) => [strand, createOutputExample(strand)])
    );

    function generateModificationSectionHeader() {
      return ui.divH([
        ui.div([ui.divText('#')], {style: {width: '20px'}}),
        ui.block75([ui.divText('Modification')]),
        ui.div([ui.divText('PTO')]),
      ]);
    }

    function generateStrandModificationPanel(strand: StrandType) {
      return ui.block([
        ui.h1(`${STRAND_LABEL[strand]}`),
        generateModificationSectionHeader(),
        modificationControls[strand],
      ], {style: {paddingTop: '12px'}});
    }

    const modificationSection = Object.fromEntries(
      STRANDS.map((strand) => [strand, generateStrandModificationPanel(strand)])
    );

    function styleExampleInputField(exampleInput: StringInput) {
      Object.assign(exampleInput.input.style, {
        resize: 'none',
        minWidth: 'none',
        flexGrow: '1',
      });
    }

    function appendOptionsToOutputInput(output: StringInput) {
      const options = ui.div([
        ui.button(ui.iconFA('copy', () => {}), () => {
          navigator.clipboard.writeText(output.value!).then(() =>
            grok.shell.info('Sequence was copied to clipboard'));
        }),
      ], 'ui-input-options');

      options.style.height = 'inherit';
      output.root.append(options);
    }

    STRANDS.forEach((strand) => {
      styleExampleInputField(inputExample[strand]);
      styleExampleInputField(outputExample[strand]);
      appendOptionsToOutputInput(outputExample[strand]);
    });



    // const inputIdColumnDiv = ui.div([]);
    const svgDisplayDiv = ui.div([]);
    const asExampleDiv = ui.div([], 'ui-form ui-form-wide');
    const loadPatternDiv = ui.div([]);
    const asModificationDiv = ui.form([]);

    function updateListOfModificationsWithNumericLabels(value: boolean) {
      if (value) {
        if (!modificationsWithNumericLabels.includes(defaultNucleotideBase)) {
          modificationsWithNumericLabels.push(defaultNucleotideBase);
        }
      } else {
        const index = modificationsWithNumericLabels.indexOf(defaultNucleotideBase, 0);
        if (index > -1) {
          modificationsWithNumericLabels.splice(index, 1);
        }
      }
    }

    const numericLabelTogglesContainer = ui.divH([
      ui.boolInput(defaultNucleotideBase, true, (value: boolean) => {
        updateListOfModificationsWithNumericLabels(value);
        refreshSvgDisplay();
        refreshOutputExamples();
      }).root,
    ]);

    const asLengthDiv = ui.div([strandLengthInputs[STRAND.ANTISENSE].root]);

    function getTableInput(tableList: DG.DataFrame[]): DG.InputBase<DG.DataFrame | null> {
      function updateStrandColumns(table: DG.DataFrame) {
        const columnNames = table.columns.names();
        STRANDS.forEach((strand) => {
          const defaultColumn = columnNames[0];
          validateStrandColumn(defaultColumn, strand);
          selectedStrandColumn[strand] = defaultColumn;
          const input = ui.choiceInput(`${STRAND_LABEL[strand]} column`, defaultColumn, columnNames, (colName: string) => {
            validateStrandColumn(colName, strand);
            selectedStrandColumn[strand] = colName;
          });
          $(strandColumnInput[strand].root).replaceWith(input.root);
        });
      }

      function updateIdColumnInput(columnNames: string[]) {
        selectedIdColumn = columnNames[0];
        const idInput = ui.choiceInput('ID column', selectedIdColumn, columnNames, (colName: string) => {
          validateIdsColumn(colName);
          selectedIdColumn = colName;
        });
        $(idColumnSelector.root).replaceWith(idInput.root);
      }

      function addTableViewIfNeeded(table: DG.DataFrame) {
        const tableName = table.name;
        if (!grok.shell.tableNames.includes(tableName)) {
          const previousView = grok.shell.v;
          grok.shell.addTableView(table);
          grok.shell.v = previousView;
        }
      }

      const tableInput = ui.tableInput('Tables', tableList[0], tableList, () => {
        const table = tableInput.value;
        if (!table) {
          console.warn('Table is null');
          return;
        }
        addTableViewIfNeeded(table);
        updateStrandColumns(table);
        updateIdColumnInput(table.columns.names());
      });

      return tableInput;
    }

    const tableInput = getTableInput([]);

    let selectedIdColumn = '';
    const idColumnSelector = ui.choiceInput('ID column', '', [], (colName: string) => {
      validateIdsColumn(colName);
      selectedIdColumn = colName;
    });
    // inputIdColumnDiv.append(inputIdColumn.root);

    refreshPatternsList();

    function toggleUiElementsBasedOnAsStrand(value: boolean) {
      const elementsToToggle = [
        modificationSection[STRAND.ANTISENSE],
        // strandColumnInputDiv[STRAND.ANTISENSE],
        strandColumnInput[STRAND.ANTISENSE].root,
        asLengthDiv,
        asModificationDiv,
        asExampleDiv,
        firstPto[STRAND.ANTISENSE].root
      ];

      elementsToToggle.forEach(element => {
        element.hidden = !value;
      });

      refreshSvgDisplay();
    }

    // Create the boolean input UI component with the refactored event handler.
    const createAsStrand = ui.boolInput('Anti sense strand', true, toggleUiElementsBasedOnAsStrand);
    createAsStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    const patternNameInput = ui.textInput('Save as', 'Pattern name', () => refreshSvgDisplay());
    patternNameInput.setTooltip('Name Of New Pattern');


    TERMINI.forEach((terminal) => {
      asModificationDiv.append(terminalModificationInputs[STRAND.ANTISENSE][terminal].root);
    })

    const patternCommentInput = ui.textInput('Comment', '', () => refreshSvgDisplay());

    function processSaveButtonClick() {
      if (patternNameInput.value !== '') {
        savePatternWithName(patternNameInput.value);
      } else {
        requestPatternNameAndSave();
      }
    }

    function savePatternWithName(patternName: string) {
      patternNameInput.value = patternName;
      savePatternInUserDataStorage().then(() => grok.shell.info('Pattern saved'));
    }

    function requestPatternNameAndSave() {
      const nameInput = ui.stringInput('Enter name', '');
      ui.dialog('Pattern Name')
        .add(nameInput.root)
        .onOK(() => savePatternWithName(nameInput.value))
        .show();
    }

    const savePatternButton = ui.bigButton('Save', processSaveButtonClick);
    patternNameInput.addOptions(savePatternButton);

    function checkForRequiredColumnSelection() {
      const condition = [true, createAsStrand.value];
      return STRANDS.some((strand, i) => condition[i] && selectedStrandColumn[strand] === '');
    }

    function checkForLengthMismatch() {
      return STRANDS.some((strand) => strandLengthInputs[strand].value !== inputExample[strand].value.length);
    }

    function adjustStrandLengthsFromTableData() {
      STRANDS.forEach((strand) => {
        strandLengthInputs[strand].value = tableInput.value!.getCol(strandColumnInput[strand].value!).getString(0).length;
      });
    }

    function addColumnsWithTranslatedSequences() {
      if (selectedIdColumn !== '') {
        addColumnWithIds(tableInput.value!.name, selectedIdColumn, getShortName(patternNameInput.value));
      }
      const condition = [true, createAsStrand.value];
      STRANDS.forEach((strand, i) => {
        if (condition[i]) {
          addColumnWithTranslatedSequences(
            tableInput.value!.name, selectedStrandColumn[strand], nucleobaseInputs[strand], ptoLinkageInputs[strand],
            terminalModificationInputs[strand][TERMINUS.FIVE_PRIME], terminalModificationInputs[strand][TERMINUS.THREE_PRIME], firstPto[strand].value!
          );
        }
      });
      grok.shell.v = grok.shell.getTableView(tableInput.value!.name);
      const columnPhrase = createAsStrand.value ? 'Columns were' : 'Column was';
      grok.shell.info(`${columnPhrase} added to table '${tableInput.value!.name}'`);
    }

    function processConvertButtonClick() {
      if (checkForRequiredColumnSelection()) {
        grok.shell.info('Please select table and columns on which to apply pattern');
        return;
      }

      if (checkForLengthMismatch()) {
        const dialog = ui.dialog('Length Mismatch');
        $(dialog.getButton('OK')).hide();
        dialog
          .add(ui.divText('Length of sequences in columns doesn\'t match entered length. Update length value?'))
          .addButton('YES', () => {
            adjustStrandLengthsFromTableData();
            dialog.close();
          })
          .show();
        return;
      }

      addColumnsWithTranslatedSequences();
      refreshOutputExamples();
    }

    const convertSequenceButton = ui.bigButton('Convert', processConvertButtonClick);


    asExampleDiv.append(inputExample[STRAND.ANTISENSE].root);
    asExampleDiv.append(outputExample[STRAND.ANTISENSE].root);

    refreshUIForNewSequenceLength();

    const svgDownloadLink = ui.link('Download', () => svgExport.saveSvgAsPng(document.getElementById('mySvg'), patternNameInput.value,
      {backgroundColor: 'white'}), 'Download pattern as PNG image', '');

    const editPatternLink = ui.link('Edit pattern', ()=>{
      ui.dialog('Edit pattern')
        .add(ui.divV([
          ui.h1('PTO'),
          ui.divH([
            isFullyPtoInput.root,
            firstPto[STRAND.SENSE].root,
            firstPto[STRAND.ANTISENSE].root,
          ], {style:{gap:'12px'}})
        ]))
        .add(ui.divH([
          modificationSection[STRAND.SENSE],
          modificationSection[STRAND.ANTISENSE],
        ], {style:{gap:'24px'}}))
        .onOK(()=>{grok.shell.info('Saved')})
        .show()
    }, 'Edit pattern', '');

    strandLengthInputs[STRAND.SENSE].addCaption('Length');

    function createLeftPanelLayout() {
      return ui.box(
        ui.div([
          ui.h1('Pattern'),
          createAsStrand.root,
          strandLengthInputs[STRAND.SENSE],
          strandLengthInputs[STRAND.ANTISENSE],
          sequenceBase.root,
          patternCommentInput.root,
          loadPatternDiv,
          patternNameInput.root,
          ui.h1('Convert'),
          tableInput.root,
          strandColumnInput[STRAND.SENSE],
          strandColumnInput[STRAND.ANTISENSE],
          idColumnSelector.root,
          ui.buttonsInput([convertSequenceButton]),
        ], 'ui-form'),
        {style: {maxWidth: '450px'}}
      );
    }

    function createRightPanelLayout() {
      return ui.panel([
        svgDisplayDiv,
        numericLabelTogglesContainer,
        generateDownloadAndEditControls(),
        generateStrandSectionDisplays(),
        ui.h1('Additional modifications'),
        ui.form([
          terminalModificationInputs[STRAND.SENSE][TERMINUS.FIVE_PRIME],
          terminalModificationInputs[STRAND.SENSE][TERMINUS.THREE_PRIME],
        ]),
        asModificationDiv,
      ], {style: {overflowX: 'scroll', padding: '12px 24px'}});
    }

    function generateDownloadAndEditControls() {
      return ui.divH([
        svgDownloadLink,
        editPatternLink
      ], {style: {gap: '12px', marginTop: '12px'}});
    }

    function generateStrandSectionDisplays() {
      return ui.divH([
        ui.divV([
          ui.h1('Sense strand'),
          inputExample[STRAND.SENSE].root,
          outputExample[STRAND.SENSE].root,
        ], 'ui-block'),
        ui.divV([
          ui.h1('Anti sense'),
          inputExample[STRAND.ANTISENSE].root,
          outputExample[STRAND.ANTISENSE].root,
        ], 'ui-block'),
      ], {style: {gap: '24px', marginTop: '24px'}});
    }

    return ui.splitH([
      createLeftPanelLayout(),
      createRightPanelLayout()
    ], {}, true);
  }
}
