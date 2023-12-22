/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsStyleMap} from '../../model/data-loading-utils/json-loader';
import {
  DEFAULT_PTO, DEFAULT_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, SS, AS, STRAND_NAME, STRANDS, TERMINAL_KEYS, TERMINAL, THREE_PRIME, FIVE_PRIME, JSON_FIELD as FIELD, StrandType, TerminalType
} from '../../model/pattern-app/const';
import {isOverhang} from '../../model/pattern-app/helpers';
import {generateExample, translateSequence, getShortName, isCurrentUserCreatedThisPattern, findDuplicates, addColumnWithIds, addColumnWithTranslatedSequences} from '../../model//pattern-app/oligo-pattern';
import {drawAxolabsPattern} from '../../model/pattern-app/draw-svg';
// todo: remove ts-ignore
//@ts-ignore
import * as svg from 'save-svg-as-png';
import $ from 'cash-dom';

const enum UPDATE_TYPE {
  PTO,
  BASIS
}

interface PatternData {
  [FIELD.SS_BASES]: string[];
  [FIELD.AS_BASES]: string[];
  [FIELD.SS_PTO]:  boolean[];
  [FIELD.AS_PTO]: boolean[];
  [FIELD.SS_3]: string;
  [FIELD.SS_5]: string;
  [FIELD.AS_3]: string;
  [FIELD.AS_5]: string;
  [FIELD.COMMENT]: string;
}

type BooleanInput = DG.InputBase<boolean | null>;
type StringInput = DG.InputBase<string | null>;

export class PatternLayoutHandler {
  get htmlDivElement() {
    function updateStrandModificationItems(strand: string) {
      clearModificationItems(strand);
      extendStrandInputArrays(strand);
      let nucleotideCounter = 0;

      for (let i = 0; i < getStrandLength(strand); i++) {
        updatePtoLinkageInputs(strand, i);
        updateBaseInputs(strand, i);

        if (!isOverhang(getBaseInputValue(strand, i))) {
          nucleotideCounter++;
        }

        const modificationItem = createModificationItem(strand, nucleotideCounter, i);
        modificationItems[strand].append(modificationItem);
      }
    }

    function clearModificationItems(strand: string) {
      modificationItems[strand].innerHTML = '';
    }

    function extendStrandInputArrays(strand: string) {
      const extensionLength = maxStrandLength[strand] - baseInputs[strand].length;
      ptoLinkageInputs[strand] = ptoLinkageInputs[strand].concat(Array(extensionLength).fill(fullyPto));
      baseInputs[strand] = baseInputs[strand].concat(Array(extensionLength).fill(sequenceBase));
    }

    function getStrandLength(strand: string) {
      return strandLengthInput[strand].value!;
    }

    function updatePtoLinkageInputs(strand: string, index: number) {
      ptoLinkageInputs[strand][index] = ui.boolInput('', ptoLinkageInputs[strand][index].value!, () => {
        refreshSvgDisplay();
        refreshOutputExamples();
      });
    }

    function updateBaseInputs(strand: string, index: number) {
      baseInputs[strand][index] = ui.choiceInput('', getBaseInputValue(strand, index), baseChoices, (v: string) => {
        handleBaseInputChange(v);
        updateStrandModificationItems(AS);
        refreshSvgDisplay();
        refreshOutputExamples();
      });

      $(baseInputs[strand][index].root).addClass('st-pattern-choice-input');
    }

    function getBaseInputValue(strand: string, index: number) {
      return baseInputs[strand][index].value!;
    }

    function handleBaseInputChange(value: string) {
      if (!enumerateModifications.includes(value)) {
        toggleModificationInList(value);
      }
    }

    function updateModificationsList(value: string, shouldAdd: boolean): void {
      const index = enumerateModifications.indexOf(value);
      if (shouldAdd && index === -1) {
        enumerateModifications.push(value);
      } else if (!shouldAdd && index > -1) {
        enumerateModifications.splice(index, 1);
      }
    }

    function appendToNumberedModificationsContainer(
      value: string,
      onChangeCallback: (value: string, shouldAdd: boolean) => void
    ): void {
      const boolInput = ui.boolInput(value, true, (boolV: boolean) => onChangeCallback(value, boolV));
      numberedModificationsListDiv.append(
        ui.divText('', {style: {width: '25px'}}),
        boolInput.root,
      );
    }

    function toggleModificationInList(value: string): void {
      const shouldAdd = !enumerateModifications.includes(value);
      updateModificationsList(value, shouldAdd);
      if (shouldAdd) {
        appendToNumberedModificationsContainer(value, (val, shouldAdd) => {
          updateModificationsList(val, shouldAdd);
          refreshSvgDisplay();
        });
      }
    }

    function createModificationItem(strand: string, nucleotideCounter: number, index: number) {
      const labelText = isOverhang(getBaseInputValue(strand, index)) ? '' : String(nucleotideCounter);
      const labelUI = ui.div([ui.label(labelText)], {style: {width: '20px'}});
      const baseInputUI = ui.block75([baseInputs[strand][index].root]);
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
      return Object.values(strandLengthInput).every(input => input.value! < MAX_SEQUENCE_LENGTH);
    }

    function extendStrandsToNewLength() {
      STRANDS.forEach(strand => {
        if (strandLengthInput[strand].value! > maxStrandLength[strand]) {
          maxStrandLength[strand] = strandLengthInput[strand].value!;
        }
        updateStrandModificationItems(strand);
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
      Object.values(strandLengthInput).forEach(input => input.value = MAX_SEQUENCE_LENGTH);
    }

    function updateInputValues(type: UPDATE_TYPE, newValue: boolean | string): void {
      const inputs = type === UPDATE_TYPE.PTO ? ptoLinkageInputs : baseInputs;

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
          inputExample[strand].value = generateExample(strandLengthInput[strand].value!, sequenceBase.value!);
        }
      });
    }

    function refreshOutputExamples(): void {
      const conditions = [true, createAsStrand.value];
      STRANDS.forEach((strand, idx) => {
        if (conditions[idx]) {
          outputExample[strand].value = translateSequence(
            inputExample[strand].value,
            baseInputs[strand],
            ptoLinkageInputs[strand],
            terminalModification[strand][FIVE_PRIME],
            terminalModification[strand][THREE_PRIME],
            firstPto[strand].value!
          );
        }
      });
    }

    function getBaseInputValues(strand: string) {
      return baseInputs[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value!);
    }

    function getPtoLinkageValues(strand: string): boolean[] {
      return [firstPto[strand].value!, ...ptoLinkageInputs[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value!)];
    }

    function refreshSvgDisplay() {
      svgDiv.innerHTML = '';
      const baseInputValues = Object.fromEntries(STRANDS.map((strand) => [strand, getBaseInputValues(strand)]));
      const ptoLinkageValues = Object.fromEntries(STRANDS.map((strand) => [strand, getPtoLinkageValues(strand)]));

      svgDiv.append(
        ui.span([
          // todo: refactor the funciton, reduce # of args
          drawAxolabsPattern(
            getShortName(saveAs.value),
            createAsStrand.value!,

            baseInputValues[SS],
            baseInputValues[AS],

            ptoLinkageValues[SS],
            ptoLinkageValues[AS],

            terminalModification[SS][THREE_PRIME].value!,
            terminalModification[SS][FIVE_PRIME].value!,

            terminalModification[AS][THREE_PRIME].value!,
            terminalModification[AS][FIVE_PRIME].value!,

            comment.value,
            enumerateModifications,
          ),
        ]),
      );
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

    async function fetchPatternFromStorage(newName: string): Promise<PatternData> {
      const entities = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);
      return JSON.parse(entities[newName]);
    }

    async function applyPatternDataToUI(newName: string, obj: PatternData): Promise<void> {
      sequenceBase.value = detectMostFrequentBase([...obj[FIELD.AS_BASES], ...obj[FIELD.SS_BASES]]);
      createAsStrand.value = (obj[FIELD.AS_BASES].length > 0);
      saveAs.value = newName;

      refreshBaseInputsFromPattern(obj);
      refreshPtoLinkagesFromPattern(obj);
      adjustStrandLengthsFromPattern(obj);
      refreshTerminalModificationsFromPattern(obj);

      comment.value = obj[FIELD.COMMENT];
    }

    async function fetchAndUpdatePatternInUI(newName: string): Promise<void> {
      const pi = DG.TaskBarProgressIndicator.create('Loading pattern...');
      try {
        const patternObj = await fetchPatternFromStorage(newName);
        await applyPatternDataToUI(newName, patternObj);
      } catch (error) {
        console.error("Error parsing pattern and updating UI: ", error);
      } finally {
        pi.close();
      }
    }

    function refreshBaseInputsFromPattern(obj: PatternData): void {
      const fields = [FIELD.SS_BASES, FIELD.AS_BASES];
      STRANDS.forEach((strand, i) => {
        baseInputs[strand] = (obj[fields[i]] as string[]).map((base: string) => ui.choiceInput('', base, baseChoices));
      });
    }

    function refreshPtoLinkagesFromPattern(obj: PatternData): void {
      const fields = [FIELD.SS_PTO, FIELD.AS_PTO];
      STRANDS.forEach((strand, i) => {
        const ptoValues = obj[fields[i]] as boolean[];
        firstPto[strand].value = ptoValues[0];
        ptoLinkageInputs[strand] = ptoValues.slice(1).map((value: boolean) => ui.boolInput('', value));
      });
    }

    function adjustStrandLengthsFromPattern(obj: PatternData): void {
      const fields = [FIELD.SS_BASES, FIELD.AS_BASES];
      STRANDS.forEach((strand, i) => {
        strandLengthInput[strand].value = obj[fields[i]].length;
      });
    }

    function refreshTerminalModificationsFromPattern(obj: PatternData): void {
      const field = [[FIELD.SS_3, FIELD.SS_5], [FIELD.AS_3, FIELD.AS_5]];
      STRANDS.forEach((strand, i) => {
        TERMINAL_KEYS.forEach((terminal, j) => {
          terminalModification[strand][terminal].value = obj[field[i][j]] as string;
        });
      });
    }

    function verifyUniformColumnLengths(colName: string): boolean {
      const col = tableInput.value!.getCol(colName);
      const areLengthsUniform = col.toList().every((value, index, array) =>
        index === 0 || value.length === array[index - 1].length || value.length === 0);

      if (!areLengthsUniform) {
        displayLengthMismatchDialog(colName);
      } else if (col.get(0).length !== strandLengthInput[SS].value) {
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

    async function savePatternAndNotify(): Promise<void> {
      const currUserName = await fetchCurrentUserName();
      saveAs.value = generatePatternSaveName(currUserName);

      const patternData = assemblePatternData();
      await grok.dapi.userDataStorage.postValue(USER_STORAGE_KEY, saveAs.value, JSON.stringify(patternData), false);

      grok.shell.info(`Pattern '${saveAs.value}' was successfully uploaded!`);
    }

    function generatePatternSaveName(currUserName: string): string {
      return saveAs.stringValue.includes('(created by ') ? 
        `${getShortName(saveAs.value)}${currUserName}` : 
        `${saveAs.stringValue}${currUserName}`;
    }

    function assemblePatternData(): PatternData {
      const createBasesArray = (strand: string) => baseInputs[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value) as string[];
      const createPtoArray = (strand: string) => [firstPto[strand].value, ...ptoLinkageInputs[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value)] as boolean[];

      return {
        [FIELD.SS_BASES]: createBasesArray(SS),
        [FIELD.AS_BASES]: createBasesArray(AS),
        [FIELD.SS_PTO]: createPtoArray(SS),
        [FIELD.AS_PTO]: createPtoArray(AS),
        [FIELD.SS_3]: terminalModification[SS][THREE_PRIME].value!,
        [FIELD.SS_5]: terminalModification[SS][FIVE_PRIME].value!,
        [FIELD.AS_3]: terminalModification[AS][THREE_PRIME].value!,
        [FIELD.AS_5]: terminalModification[AS][FIVE_PRIME].value!,
        [FIELD.COMMENT]: comment.value,
      };
    }

    async function sortPatternsByOwnership(patternsData: PatternData): Promise<{ lstMy: string[], lstOthers: string[] }> {
      const lstMy: string[] = [];
      const lstOthers: string[] = [];

      for (const ent of Object.keys(patternsData)) {
        if (await isCurrentUserCreatedThisPattern(ent))
          lstOthers.push(ent);
        else
          lstMy.push(ent);
      }

      return { lstMy, lstOthers };
    }

    function initializeLoadPatternInterface(loadPattern: StringInput, loadPatternDiv: HTMLElement, patternListChoiceInput: StringInput) {
      clearAndRepopulateUIElements(loadPatternDiv, loadPattern, patternListChoiceInput);
      applyStyleToLoadPatternInput(loadPattern.input);
      loadPattern.setTooltip('Apply Existing Pattern');
      const deleteButton = generateDeletePatternButton(loadPattern);
      loadPattern.root.append(deleteButton);
    }

    function clearAndRepopulateUIElements(loadPatternDiv: HTMLElement, loadPattern: StringInput, patternListChoiceInput: StringInput) {
      loadPatternDiv.innerHTML = '';
      loadPattern.root.append(patternListChoiceInput.input, loadPattern.input);
      loadPatternDiv.append(loadPattern.root);
    }

    function applyStyleToLoadPatternInput(loadPatternInput: HTMLElement) {
      loadPatternInput.style.maxWidth = '120px';
      loadPatternInput.style.marginLeft = '12px';
    }

    function generateDeletePatternButton(loadPattern: StringInput): HTMLElement {
      return ui.div([
        ui.button(ui.iconFA('trash-alt'), async () => {
          if (!loadPattern.value) {
            grok.shell.warning('Choose pattern to delete');
          } else if (await isCurrentUserCreatedThisPattern(saveAs.value)) {
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
      const { lstMy, lstOthers } = await sortPatternsByOwnership(patternsData);

      const currentUserName = (await grok.dapi.users.current()).friendlyName;
      const otherUsers = 'Other users';

      let loadPattern = ui.choiceInput('Load pattern', '', lstMy, (v: string) => fetchAndUpdatePatternInUI(v));

      const patternListChoiceInput = ui.choiceInput(
        '', currentUserName, [currentUserName, otherUsers],
        (v: string) => patternListChoiceInputOnChange(v)
      );

      patternListChoiceInput.input.style.maxWidth = '142px';
      initializeLoadPatternInterface(loadPattern, loadPatternDiv, patternListChoiceInput);

      function patternListChoiceInputOnChange(v: string) {
        const currentList = v === currentUserName ? lstMy : lstOthers;
        loadPattern = ui.choiceInput('Load pattern', '', currentList, (v: string) => fetchAndUpdatePatternInUI(v));
        initializeLoadPatternInterface(loadPattern, loadPatternDiv, patternListChoiceInput);
      }
    }

    async function checkIfPatternExistsInStorage(patternData: PatternData, patternName: string) {
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
          await savePatternAndNotify();
          dialog.close();
        })
        .show();
    }

    async function deletePatternFromStorage(patternName: string) {
      return grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, patternName, false);
    }

    async function savePatternInUserDataStorage(): Promise<void> {
      const patternData = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);

      const patternName = saveAs.value;

      if (await checkIfPatternExistsInStorage(patternData, patternName)) {
        await displayPatternReplaceConfirmationDialog(patternName);
      } else {
        await savePatternAndNotify();
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
      const currentStrandLength = strandLengthInput[strand].value;
      if (sequenceLength !== currentStrandLength) {
        strandLengthInput[strand].value = sequenceLength;
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


    const baseChoices: string[] = Object.keys(axolabsStyleMap);
    const defaultBase: string = baseChoices[0];
    const enumerateModifications = [defaultBase];
    const sequenceBase = ui.choiceInput('Sequence basis', defaultBase, baseChoices, (v: string) => {
      // updateBases(v);
      updateInputValues(UPDATE_TYPE.BASIS, v);
      refreshOutputExamples();
    });
    
    function createFullyPtoInput(): BooleanInput {
      const fullyPto = ui.boolInput('Fully PTO', DEFAULT_PTO, handleFullPtoInputChange);
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

    function handleFullPtoInputChange(v: boolean): void {
      STRANDS.forEach((s) => { firstPto[s].value = v; });
      // updatePto(v);
      updateInputValues(UPDATE_TYPE.PTO, v);
      refreshOutputExamples();
    }

    const fullyPto = createFullyPtoInput();

    const maxStrandLength = Object.fromEntries(STRANDS.map(
      (strand) => [strand, DEFAULT_SEQUENCE_LENGTH]
    ));
    // todo: remove vague legacy 'items' from name
    const modificationItems = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.div([])]
    ));
    const ptoLinkageInputs = Object.fromEntries(STRANDS.map(
      (strand) => [strand, Array<BooleanInput>(DEFAULT_SEQUENCE_LENGTH)
        .fill(ui.boolInput('', DEFAULT_PTO))]
    ));
    const baseInputs = Object.fromEntries(STRANDS.map(
      (strand) => {
        const choiceInputs = Array<StringInput>(DEFAULT_SEQUENCE_LENGTH)
          .fill(ui.choiceInput('', defaultBase, baseChoices));
        return [strand, choiceInputs];
      }
    ));
    const strandLengthInput = Object.fromEntries(STRANDS.map(
      (strand) => {
        const input = ui.intInput(`${STRAND_NAME[strand]} length`, DEFAULT_SEQUENCE_LENGTH, () => refreshUIForNewSequenceLength());
        input.setTooltip(`Length of ${STRAND_NAME[strand].toLowerCase()}, including overhangs`);
        return [strand, input];
      }));
    const strandVar = Object.fromEntries(STRANDS.map((strand) => [strand, '']));
    const inputExample = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.textInput(
        ``, generateExample(strandLengthInput[strand].value!, sequenceBase.value!))
      ]));

    const strandColumnInput = Object.fromEntries(STRANDS.map((strand) => {
      const input = ui.choiceInput(`${STRAND_NAME[strand]} column`, '', [], (colName: string) => {
        validateStrandColumn(colName, strand);
        strandVar[strand] = colName;
      });
      return [strand, input];
    }));

    function createFirstPtoInputs(fullyPto: BooleanInput): Record<string,BooleanInput> {
      return Object.fromEntries(STRANDS.map((strand) => [strand, createStrandPtoInput(strand, fullyPto)]));
    }

    function createStrandPtoInput(strand: StrandType, fullyPto: BooleanInput): BooleanInput {
      const input = ui.boolInput(`First ${strand} PTO`, fullyPto.value!, refreshSvgDisplay);
      input.setTooltip(`ps linkage before first nucleotide of ${STRAND_NAME[strand].toLowerCase()}`);
      configureInputLabel(input.captionLabel);
      return input;
    }

    function configureInputLabel(label: HTMLElement): void {
      label.classList.add('ui-label-right');
      Object.assign(label.style, {
        textAlign: 'left',
        maxWidth: '100px',
        minWidth: '40px',
        width: 'auto'
      });
    }

    const firstPto = createFirstPtoInputs(fullyPto);

    function createTerminalModificationInputs() {
      return Object.fromEntries(STRANDS.map(strand => [strand, createStrandInputs(strand)]));
    }

    function createStrandInputs(strand: StrandType) {
      return Object.fromEntries(TERMINAL_KEYS.map((key: TerminalType) => [key, createModificationInput(strand, key)]));
    }

    function createModificationInput(strand: StrandType, key: TerminalType): StringInput {
      const label = `${strand} ${TERMINAL[key]}\' Modification`;
      const tooltip = `Additional ${strand} ${TERMINAL[key]}\' Modification`;
      const input = ui.stringInput(label, '', modificationInputChangeCallback);

      input.setTooltip(tooltip);
      return input;
    }

    function modificationInputChangeCallback(): void {
      refreshSvgDisplay();
      refreshOutputExamples();
    }

    const terminalModification = createTerminalModificationInputs();

    function createOutputExample(strand: StrandType) {
      const translatedSequence = translateSequence(
        inputExample[strand].value,
        baseInputs[strand],
        ptoLinkageInputs[strand],
        terminalModification[strand][THREE_PRIME],
        terminalModification[strand][FIVE_PRIME],
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

    function createModificationHeader() {
      return ui.divH([
        ui.div([ui.divText('#')], {style: {width: '20px'}}),
        ui.block75([ui.divText('Modification')]),
        ui.div([ui.divText('PTO')]),
      ]);
    }

    function createModificationPanel(strand: StrandType) {
      return ui.block([
        ui.h1(`${STRAND_NAME[strand]}`),
        createModificationHeader(),
        modificationItems[strand],
      ], {style: {paddingTop: '12px'}});
    }

    const modificationSection = Object.fromEntries(
      STRANDS.map((strand) => [strand, createModificationPanel(strand)])
    );

    function applyExampleStyles(example: StringInput) {
      Object.assign(example.input.style, {
        resize: 'none',
        minWidth: 'none',
        flexGrow: '1',
      });
    }

    function appendOptionsToOutput(output: StringInput) {
      const options = ui.div([
        ui.button(ui.iconFA('copy', () => {}), () => {
          navigator.clipboard.writeText(output.value!).then(() =>
            grok.shell.info('Sequence was copied to clipboard'));
        }),
      ], 'ui-input-options');

      options.style.height = 'inherit';
      output.root.append(options);
    }

    STRANDS.forEach((s) => {
      applyExampleStyles(inputExample[s]);
      applyExampleStyles(outputExample[s]);
      appendOptionsToOutput(outputExample[s]);
    });



    // const inputIdColumnDiv = ui.div([]);
    const svgDiv = ui.div([]);
    const asExampleDiv = ui.div([], 'ui-form ui-form-wide');
    const loadPatternDiv = ui.div([]);
    const asModificationDiv = ui.form([]);

    function updateEnumerateModificationsList(value: boolean) {
      if (value) {
        if (!enumerateModifications.includes(defaultBase)) {
          enumerateModifications.push(defaultBase);
        }
      } else {
        const index = enumerateModifications.indexOf(defaultBase, 0);
        if (index > -1) {
          enumerateModifications.splice(index, 1);
        }
      }
    }

    const numberedModificationsListDiv = ui.divH([
      ui.boolInput(defaultBase, true, (v: boolean) => {
        updateEnumerateModificationsList(v);
        refreshSvgDisplay();
        refreshOutputExamples();
      }).root,
    ]);

    const asLengthDiv = ui.div([strandLengthInput[AS].root]);

    function getTableInput(tableList: DG.DataFrame[]): DG.InputBase<DG.DataFrame | null> {
      function updateStrandColumns(table: DG.DataFrame) {
        const columnNames = table.columns.names();
        STRANDS.forEach((strand) => {
          const defaultColumn = columnNames[0];
          validateStrandColumn(defaultColumn, strand);
          strandVar[strand] = defaultColumn;
          const input = ui.choiceInput(`${STRAND_NAME[strand]} column`, defaultColumn, columnNames, (colName: string) => {
            validateStrandColumn(colName, strand);
            strandVar[strand] = colName;
          });
          $(strandColumnInput[strand].root).replaceWith(input.root);
        });
      }

      function updateIdColumnInput(columnNames: string[]) {
        idVar = columnNames[0];
        const idInput = ui.choiceInput('ID column', idVar, columnNames, (colName: string) => {
          validateIdsColumn(colName);
          idVar = colName;
        });
        $(inputIdColumn.root).replaceWith(idInput.root);
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

    // todo: unify with strandVar
    let idVar = '';
    const inputIdColumn = ui.choiceInput('ID column', '', [], (colName: string) => {
      validateIdsColumn(colName);
      idVar = colName;
    });
    // inputIdColumnDiv.append(inputIdColumn.root);

    refreshPatternsList();

    function toggleUiElementsBasedOnAsStrand(value: boolean) {
      const elementsToToggle = [
        modificationSection[AS],
        // strandColumnInputDiv[AS],
        strandColumnInput[AS].root,
        asLengthDiv,
        asModificationDiv,
        asExampleDiv,
        firstPto[AS].root
      ];

      elementsToToggle.forEach(element => {
        element.hidden = !value;
      });

      refreshSvgDisplay();
    }

    // Create the boolean input UI component with the refactored event handler.
    const createAsStrand = ui.boolInput('Anti sense strand', true, toggleUiElementsBasedOnAsStrand);
    createAsStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    const saveAs = ui.textInput('Save as', 'Pattern name', () => refreshSvgDisplay());
    saveAs.setTooltip('Name Of New Pattern');


    TERMINAL_KEYS.forEach((terminal) => {
      asModificationDiv.append(terminalModification[AS][terminal].root);
    })

    const comment = ui.textInput('Comment', '', () => refreshSvgDisplay());

    function handleSave() {
      if (saveAs.value !== '') {
        savePatternWithName(saveAs.value);
      } else {
        promptForPatternNameAndSave();
      }
    }

    function savePatternWithName(patternName: string) {
      saveAs.value = patternName;
      savePatternInUserDataStorage().then(() => grok.shell.info('Pattern saved'));
    }

    function promptForPatternNameAndSave() {
      const nameInput = ui.stringInput('Enter name', '');
      ui.dialog('Pattern Name')
        .add(nameInput.root)
        .onOK(() => savePatternWithName(nameInput.value))
        .show();
    }

    const savePatternButton = ui.bigButton('Save', handleSave);
    saveAs.addOptions(savePatternButton);

    function areRequiredColumnsSelected() {
      const condition = [true, createAsStrand.value];
      return STRANDS.some((s, i) => condition[i] && strandVar[s] === '');
    }

    function isLengthMismatchPresent() {
      return STRANDS.some((s) => strandLengthInput[s].value !== inputExample[s].value.length);
    }

    function updateStrandLengthsBasedOnTable() {
      STRANDS.forEach((s) => {
        strandLengthInput[s].value = tableInput.value!.getCol(strandColumnInput[s].value!).getString(0).length;
      });
    }

    function addTranslatedSequenceColumns() {
      if (idVar !== '') {
        addColumnWithIds(tableInput.value!.name, idVar, getShortName(saveAs.value));
      }
      const condition = [true, createAsStrand.value];
      STRANDS.forEach((strand, i) => {
        if (condition[i]) {
          addColumnWithTranslatedSequences(
            tableInput.value!.name, strandVar[strand], baseInputs[strand], ptoLinkageInputs[strand],
            terminalModification[strand][FIVE_PRIME], terminalModification[strand][THREE_PRIME], firstPto[strand].value!
          );
        }
      });
      grok.shell.v = grok.shell.getTableView(tableInput.value!.name);
      const columnPhrase = createAsStrand.value ? 'Columns were' : 'Column was';
      grok.shell.info(`${columnPhrase} added to table '${tableInput.value!.name}'`);
    }

    function handleConvert() {
      if (areRequiredColumnsSelected()) {
        grok.shell.info('Please select table and columns on which to apply pattern');
        return;
      }

      if (isLengthMismatchPresent()) {
        const dialog = ui.dialog('Length Mismatch');
        $(dialog.getButton('OK')).hide();
        dialog
          .add(ui.divText('Length of sequences in columns doesn\'t match entered length. Update length value?'))
          .addButton('YES', () => {
            updateStrandLengthsBasedOnTable();
            dialog.close();
          })
          .show();
        return;
      }

      addTranslatedSequenceColumns();
      refreshOutputExamples();
    }

    const convertSequenceButton = ui.bigButton('Convert', handleConvert);


    asExampleDiv.append(inputExample[AS].root);
    asExampleDiv.append(outputExample[AS].root);

    refreshUIForNewSequenceLength();

    const inputsSection = ui.block50([
      ui.h1('Convert options'),
      tableInput.root,
      strandColumnInput[SS].root,
      strandColumnInput[AS].root,
      inputIdColumn.root,
      ui.buttonsInput([
        convertSequenceButton,
      ]),
    ]);
    inputsSection.classList.add('ui-form');

    const downloadButton = ui.link('Download', () => svg.saveSvgAsPng(document.getElementById('mySvg'), saveAs.value,
      {backgroundColor: 'white'}), 'Download pattern as PNG image', '');

    const editPattern = ui.link('Edit pattern', ()=>{
      ui.dialog('Edit pattern')
        .add(ui.divV([
          ui.h1('PTO'),
          ui.divH([
            fullyPto.root,
            firstPto[SS].root,
            firstPto[AS].root,
          ], {style:{gap:'12px'}})
        ]))
        .add(ui.divH([
          modificationSection[SS],
          modificationSection[AS],
        ], {style:{gap:'24px'}}))
        .onOK(()=>{grok.shell.info('Saved')})
        .show()
    }, 'Edit pattern', '');  

    strandLengthInput[SS].addCaption('Length');

    function createLeftPanel() {
      return ui.box(
        ui.div([
          ui.h1('Pattern'),
          createAsStrand.root,
          strandLengthInput[SS],
          strandLengthInput[AS],
          sequenceBase.root,
          comment.root,
          loadPatternDiv,
          saveAs.root,
          ui.h1('Convert'),
          tableInput.root,
          strandColumnInput[SS],
          strandColumnInput[AS],
          inputIdColumn.root,
          ui.buttonsInput([convertSequenceButton]),
        ], 'ui-form'),
        {style: {maxWidth: '450px'}}
      );
    }

    function createRightPanel() {
      return ui.panel([
        svgDiv,
        numberedModificationsListDiv,
        createDownloadEditButtons(),
        createStrandSections(),
        ui.h1('Additional modifications'),
        ui.form([
          terminalModification[SS][FIVE_PRIME],
          terminalModification[SS][THREE_PRIME],
        ]),
        asModificationDiv,
      ], {style: {overflowX: 'scroll', padding: '12px 24px'}});
    }

    function createDownloadEditButtons() {
      return ui.divH([
        downloadButton,
        editPattern
      ], {style: {gap: '12px', marginTop: '12px'}});
    }

    function createStrandSections() {
      return ui.divH([
        ui.divV([
          ui.h1('Sense strand'),
          inputExample[SS].root,
          outputExample[SS].root,
        ], 'ui-block'),
        ui.divV([
          ui.h1('Anti sense'),
          inputExample[AS].root,
          outputExample[AS].root,
        ], 'ui-block'),
      ], {style: {gap: '24px', marginTop: '24px'}});
    }

    return ui.splitH([
      createLeftPanel(),
      createRightPanel()
    ], {}, true);
  }
}
