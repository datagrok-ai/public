/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsStyleMap} from '../../model/data-loading-utils/json-loader';
import {
  DEFAULT_PTO, DEFAULT_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, SS, AS, STRAND_NAME, STRANDS, TERMINAL, TERMINAL_KEYS, THREE_PRIME, FIVE_PRIME, JSON_FIELD as FIELD
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
    function updateModification(strand: string) {
      clearModificationItems(strand);
      extendStrandArrays(strand);
      let nucleotideCounter = 0;

      for (let i = 0; i < getStrandLength(strand); i++) {
        updatePtoLinkage(strand, i);
        updateBaseInputObject(strand, i);

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

    function extendStrandArrays(strand: string) {
      const extensionLength = maxStrandLength[strand] - baseInputsObject[strand].length;
      ptoLinkages[strand] = ptoLinkages[strand].concat(Array(extensionLength).fill(fullyPto));
      baseInputsObject[strand] = baseInputsObject[strand].concat(Array(extensionLength).fill(sequenceBase));
    }

    function getStrandLength(strand: string) {
      return strandLengthInput[strand].value!;
    }

    function updatePtoLinkage(strand: string, index: number) {
      ptoLinkages[strand][index] = ui.boolInput('', ptoLinkages[strand][index].value!, () => {
        updateSvgScheme();
        updateOutputExamples();
      });
    }

    function updateBaseInputObject(strand: string, index: number) {
      baseInputsObject[strand][index] = ui.choiceInput('', getBaseInputValue(strand, index), baseChoices, (v: string) => {
        handleBaseInputChange(v);
        updateModification(AS);
        updateSvgScheme();
        updateOutputExamples();
      });

      $(baseInputsObject[strand][index].root).addClass('st-pattern-choice-input');
    }

    function getBaseInputValue(strand: string, index: number) {
      return baseInputsObject[strand][index].value!;
    }

    function handleBaseInputChange(value: string) {
      if (!enumerateModifications.includes(value)) {
        updateEnumerateModifications(value);
      }
    }

    function updateEnumerateModifications(value: string) {
      if (!enumerateModifications.includes(value)) {
        enumerateModifications.push(value);
        appendModificationToDiv(value);
      } else {
        const index = enumerateModifications.indexOf(value);
        if (index > -1)
          enumerateModifications.splice(index, 1);
      }
    }

    function appendModificationToDiv(value: string) {
      isEnumerateModificationsDiv.append(
        ui.divText('', {style: {width: '25px'}}),
        ui.boolInput(value, true, (boolV: boolean) => {
          if (boolV && !enumerateModifications.includes(value)) {
            enumerateModifications.push(value);
          } else if (!boolV) {
            const index = enumerateModifications.indexOf(value);
            if (index > -1) {
              enumerateModifications.splice(index, 1);
            }
          }
          updateSvgScheme();
        }).root,
      );
    }

    function createModificationItem(strand: string, nucleotideCounter: number, index: number) {
      const labelText = isOverhang(getBaseInputValue(strand, index)) ? '' : String(nucleotideCounter);
      const labelUI = ui.div([ui.label(labelText)], {style: {width: '20px'}});
      const baseInputUI = ui.block75([baseInputsObject[strand][index].root]);
      const ptoLinkageUI = ui.div([ptoLinkages[strand][index]]);

      return ui.divH([labelUI, baseInputUI, ptoLinkageUI], {style: {alignItems: 'center'}});
    }

    function updateUiForNewSequenceLength() {
      if (isSequenceLengthWithinRange()) {
        updateStrandsForNewLength();
        updateUiComponents();
      } else {
        showOutOfRangeDialog();
      }
    }

    function isSequenceLengthWithinRange() {
      return Object.values(strandLengthInput).every(input => input.value! < MAX_SEQUENCE_LENGTH);
    }

    function updateStrandsForNewLength() {
      STRANDS.forEach(strand => {
        if (strandLengthInput[strand].value! > maxStrandLength[strand]) {
          maxStrandLength[strand] = strandLengthInput[strand].value!;
        }
        updateModification(strand);
      });
    }

    function updateUiComponents() {
      updateSvgScheme();
      updateInputExamples();
      updateOutputExamples();
    }

    function showOutOfRangeDialog() {
      ui.dialog('Out of range')
        .add(ui.divText(`Sequence length should be less than ${MAX_SEQUENCE_LENGTH} due to UI constraints.`))
        .onOK(() => resetStrandLengthInputs())
        .onCancel(() => resetStrandLengthInputs())
        .showModal(false);
    }

    function resetStrandLengthInputs() {
      Object.values(strandLengthInput).forEach(input => input.value = MAX_SEQUENCE_LENGTH);
    }

    function updateValues(type: UPDATE_TYPE, newValue: boolean | string): void {
      const targetObject = type === UPDATE_TYPE.PTO ? ptoLinkages : baseInputsObject;

      for (let i = 0; i < STRANDS.length; i++) {
        const strand = STRANDS[i];
        // WARNING: replacing this with for (const ...) or .forEach() leads to a bug: some values are not updated !!!
        for (let j = 0; j < targetObject[strand].length; j++) {
          const item = targetObject[strand][j];
          item.value = newValue;
        }
      }
      updateSvgScheme();
    }

    function updateInputExamples(): void {
      STRANDS.forEach((strand) => {
        if (strandColumnInput[strand].value === '') {
          inputExample[strand].value = generateExample(strandLengthInput[strand].value!, sequenceBase.value!);
        }
      });
    }

    function updateOutputExamples(): void {
      const conditions = [true, createAsStrand.value];
      STRANDS.forEach((strand, idx) => {
        if (conditions[idx]) {
          outputExample[strand].value = translateSequence(
            inputExample[strand].value,
            baseInputsObject[strand],
            ptoLinkages[strand],
            terminalModification[strand][FIVE_PRIME],
            terminalModification[strand][THREE_PRIME],
            firstPto[strand].value!
          );
        }
      });
    }

    function getBaseInputValues(strand: string) {
      return baseInputsObject[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value!);
    }

    function getPtoLinkageValues(strand: string): boolean[] {
      return [firstPto[strand].value!, ...ptoLinkages[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value!)];
    }

    function updateSvgScheme() {
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

            terminalModification[SS][THREE_PRIME].value,
            terminalModification[SS][FIVE_PRIME].value,

            terminalModification[AS][THREE_PRIME].value,
            terminalModification[AS][FIVE_PRIME].value,

            comment.value,
            enumerateModifications,
          ),
        ]),
      );
    }

    // todo: rename
    function detectDefaultBasis(array: string[]): string {
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

    async function parsePattern(newName: string): Promise<PatternData> {
      const entities = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);
      return JSON.parse(entities[newName]);
    }

    async function updateUiWithPattern(newName: string, obj: PatternData): Promise<void> {
      sequenceBase.value = detectDefaultBasis([...obj[FIELD.AS_BASES], ...obj[FIELD.SS_BASES]]);
      createAsStrand.value = (obj[FIELD.AS_BASES].length > 0);
      saveAs.value = newName;

      updateBaseInputs(obj);
      updatePtoLinkages(obj);
      updateStrandLengths(obj);
      updateTerminalModifications(obj);

      comment.value = obj[FIELD.COMMENT];
    }

    async function parsePatternAndUpdateUi(newName: string): Promise<void> {
      const pi = DG.TaskBarProgressIndicator.create('Loading pattern...');
      try {
        const patternObj = await parsePattern(newName);
        await updateUiWithPattern(newName, patternObj);
      } catch (error) {
        console.error("Error parsing pattern and updating UI: ", error);
      } finally {
        pi.close();
      }
    }

    function updateBaseInputs(obj: PatternData): void {
      const fields = [FIELD.SS_BASES, FIELD.AS_BASES];
      STRANDS.forEach((strand, i) => {
        baseInputsObject[strand] = (obj[fields[i]] as string[]).map((base: string) => ui.choiceInput('', base, baseChoices));
      });
    }

    function updatePtoLinkages(obj: PatternData): void {
      const fields = [FIELD.SS_PTO, FIELD.AS_PTO];
      STRANDS.forEach((strand, i) => {
        const ptoValues = obj[fields[i]] as boolean[];
        firstPto[strand].value = ptoValues[0];
        ptoLinkages[strand] = ptoValues.slice(1).map((value: boolean) => ui.boolInput('', value));
      });
    }

    function updateStrandLengths(obj: PatternData): void {
      const fields = [FIELD.SS_BASES, FIELD.AS_BASES];
      STRANDS.forEach((strand, i) => {
        strandLengthInput[strand].value = obj[fields[i]].length;
      });
    }

    function updateTerminalModifications(obj: PatternData): void {
      const field = [[FIELD.SS_3, FIELD.SS_5], [FIELD.AS_3, FIELD.AS_5]];
      STRANDS.forEach((strand, i) => {
        TERMINAL_KEYS.forEach((terminal, j) => {
          terminalModification[strand][terminal].value = obj[field[i][j]] as string;
        });
      });
    }

    function checkColumnLengthsUniform(colName: string): boolean {
      const col = tableInput.value!.getCol(colName);
      const areLengthsUniform = col.toList().every((value, index, array) =>
        index === 0 || value.length === array[index - 1].length || value.length === 0);

      if (!areLengthsUniform) {
        showSequencesLengthMismatchDialog(colName);
      } else if (col.get(0).length !== strandLengthInput[SS].value) {
        showLengthUpdatedDialog();
      }

      return areLengthsUniform;
    }

    function showSequencesLengthMismatchDialog(colName: string) {
      const dialog = ui.dialog('Sequences lengths mismatch');
      $(dialog.getButton('OK')).hide();

      dialog
        .add(ui.divText('The sequence length should match the number of Raw sequences in the input file'))
        .add(ui.divText('\'ADD COLUMN\' to see sequences lengths'))
        .addButton('ADD COLUMN', () => addLengthColumn(colName, dialog))
        .show();
    }

    function addLengthColumn(colName: string, dialog: any) {
      const table = tableInput.value!;
      table.columns.addNewInt('Sequences lengths in ' + colName).init((j: number) => table.getCol(colName).get(j).length);
      grok.shell.info('Column with lengths added to \'' + table.name + '\'');
      dialog.close();
      grok.shell.v = grok.shell.getTableView(table.name);
    }

    function showLengthUpdatedDialog() {
      const d = ui.dialog('Length was updated by value from imported file');
      d.add(ui.divText('Latest modifications may not take effect during translation'))
        .onOK(() => grok.shell.info('Lengths changed'))
        .show();
    }

    async function getCurrentUserName(): Promise<string> {
      const user = await grok.dapi.users.current();
      return ` (created by ${user.friendlyName})`;
    }

    async function postPatternToUserStorage(): Promise<void> {
      const currUserName = await getCurrentUserName();
      saveAs.value = createSaveAsValue(currUserName);

      const patternData = createPatternData();
      await grok.dapi.userDataStorage.postValue(USER_STORAGE_KEY, saveAs.value, JSON.stringify(patternData), false);

      grok.shell.info(`Pattern '${saveAs.value}' was successfully uploaded!`);
    }

    function createSaveAsValue(currUserName: string): string {
      return saveAs.stringValue.includes('(created by ') ? 
        `${getShortName(saveAs.value)}${currUserName}` : 
        `${saveAs.stringValue}${currUserName}`;
    }

    function createPatternData(): PatternData {
      const createBasesArray = (strand: string) => baseInputsObject[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value) as string[];
      const createPtoArray = (strand: string) => [firstPto[strand].value, ...ptoLinkages[strand].slice(0, strandLengthInput[strand].value!).map(e => e.value)] as boolean[];

      return {
        [FIELD.SS_BASES]: createBasesArray(SS),
        [FIELD.AS_BASES]: createBasesArray(AS),
        [FIELD.SS_PTO]: createPtoArray(SS),
        [FIELD.AS_PTO]: createPtoArray(AS),
        [FIELD.SS_3]: terminalModification[SS][THREE_PRIME].value,
        [FIELD.SS_5]: terminalModification[SS][FIVE_PRIME].value,
        [FIELD.AS_3]: terminalModification[AS][THREE_PRIME].value,
        [FIELD.AS_5]: terminalModification[AS][FIVE_PRIME].value,
        [FIELD.COMMENT]: comment.value,
      };
    }



    async function updatePatternsList() {
      const patternsData = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);
      const { lstMy, lstOthers } = await categorizePatterns();

      const currentUserName = (await grok.dapi.users.current()).friendlyName;
      const otherUsers = 'Other users';

      let loadPattern = ui.choiceInput('Load pattern', '', lstMy, (v: string) => parsePatternAndUpdateUi(v));
      const patternListChoiceInput = ui.choiceInput('', currentUserName, [currentUserName, otherUsers], (v: string) => {
        const currentList = v === currentUserName ? lstMy : lstOthers;
        loadPattern = ui.choiceInput('Load pattern', '', currentList, (v: string) => parsePatternAndUpdateUi(v));
        setupLoadPatternUI(loadPattern, loadPatternDiv, patternListChoiceInput);
      });

      patternListChoiceInput.input.style.maxWidth = '142px';
      setupLoadPatternUI(loadPattern, loadPatternDiv, patternListChoiceInput);

      async function categorizePatterns(): Promise<{ lstMy: string[], lstOthers: string[] }> {
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


      function setupLoadPatternUI(loadPattern: StringInput, loadPatternDiv: HTMLElement, patternListChoiceInput: StringInput) {
        clearAndAppendNewElements(loadPatternDiv, loadPattern, patternListChoiceInput);
        styleLoadPatternInput(loadPattern.input);
        loadPattern.setTooltip('Apply Existing Pattern');
        const deleteButton = createDeleteButton(loadPattern);
        loadPattern.root.append(deleteButton);
      }

      function clearAndAppendNewElements(loadPatternDiv: HTMLElement, loadPattern: StringInput, patternListChoiceInput: StringInput) {
        loadPatternDiv.innerHTML = '';
        loadPattern.root.append(patternListChoiceInput.input, loadPattern.input);
        loadPatternDiv.append(loadPattern.root);
      }

      function styleLoadPatternInput(loadPatternInput: HTMLElement) {
        loadPatternInput.style.maxWidth = '120px';
        loadPatternInput.style.marginLeft = '12px';
      }

      function createDeleteButton(loadPattern: StringInput): HTMLElement {
        return ui.div([
          ui.button(ui.iconFA('trash-alt'), async () => {
            if (!loadPattern.value) {
              grok.shell.warning('Choose pattern to delete');
            } else if (await isCurrentUserCreatedThisPattern(saveAs.value)) {
              grok.shell.warning('Cannot delete pattern, created by other user');
            } else {
              await deletePattern(loadPattern.value);
              await updatePatternsList();
            }
          })
        ], 'ui-input-options');
      }

      async function deletePattern(patternName: string) {
        await grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, patternName, false);
        grok.shell.info(`Pattern '${patternName}' deleted`);
      }
    }




    async function savePattern() {
      await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false)
        .then((entities) => {
          if (Object.keys(entities).includes(saveAs.value)) {
            const dialog = ui.dialog('Pattern already exists');
            $(dialog.getButton('OK')).hide();
            dialog
              .add(ui.divText('Pattern name \'' + saveAs.value + '\' already exists.'))
              .add(ui.divText('Replace pattern?'))
              .addButton('YES', async () => {
                await grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, saveAs.value, false)
                  .then(() => postPatternToUserStorage());
                dialog.close();
              })
              .show();
          } else
            postPatternToUserStorage();
        });
      await updatePatternsList();
    }

    function validateStrandColumn(colName: string, strand: string): void {
      const allLengthsAreTheSame: boolean = checkColumnLengthsUniform(colName);
      const firstSequence = tableInput.value!.getCol(colName).get(0);
      if (allLengthsAreTheSame && firstSequence.length !== strandLengthInput[strand].value)
      strandLengthInput[strand].value = tableInput.value!.getCol(colName).get(0).length;
      inputExample[strand].value = firstSequence;
    }

    function validateIdsColumn(colName: string) {
      const col = tableInput.value!.getCol(colName);
      if (col.type !== DG.TYPE.INT)
        grok.shell.error('Column should contain integers only');
      //@ts-ignore
      else if (col.categories.filter((e) => e !== '').length < col.toList().filter((e) => e !== '').length) {
        const duplicates = findDuplicates(col.getRawData());
        ui.dialog('Non-unique IDs')
          .add(ui.divText('Press \'OK\' to select rows with non-unique values'))
          .onOK(() => {
            const selection = tableInput.value!.selection;
            selection.init((i: number) => duplicates.indexOf(col.get(i)) > -1);
            grok.shell.v = grok.shell.getTableView(tableInput.value!.name);
            grok.shell.info('Rows are selected in table \'' + tableInput.value!.name + '\'');
          })
          .show();
      }
    }

    const baseChoices: string[] = Object.keys(axolabsStyleMap);
    const defaultBase: string = baseChoices[0];
    const enumerateModifications = [defaultBase];
    const sequenceBase = ui.choiceInput('Sequence basis', defaultBase, baseChoices, (v: string) => {
      // updateBases(v);
      updateValues(UPDATE_TYPE.BASIS, v);
      updateOutputExamples();
    });
    const fullyPto = ui.boolInput('Fully PTO', DEFAULT_PTO, (v: boolean) => {
      STRANDS.forEach((s) => { firstPto[s].value = v; })
      // updatePto(v);
      updateValues(UPDATE_TYPE.PTO, v);
      updateOutputExamples();
    });
    fullyPto.captionLabel.classList.add('ui-label-right');
    fullyPto.captionLabel.style.textAlign = 'left';
    fullyPto.captionLabel.style.maxWidth = '100px';
    fullyPto.captionLabel.style.maxWidth = '100px';
    fullyPto.captionLabel.style.minWidth = '40px';
    fullyPto.captionLabel.style.width = 'auto';

    const maxStrandLength = Object.fromEntries(STRANDS.map(
      (strand) => [strand, DEFAULT_SEQUENCE_LENGTH]
    ));
    // todo: remove vague legacy 'items' from name
    const modificationItems = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.div([])]
    ));
    const ptoLinkages = Object.fromEntries(STRANDS.map(
      (strand) => [strand, Array<BooleanInput>(DEFAULT_SEQUENCE_LENGTH)
        .fill(ui.boolInput('', DEFAULT_PTO))]
    ));
    const baseInputsObject = Object.fromEntries(STRANDS.map(
      (strand) => {
        const choiceInputs = Array<StringInput>(DEFAULT_SEQUENCE_LENGTH)
          .fill(ui.choiceInput('', defaultBase, baseChoices));
        return [strand, choiceInputs];
      }
    ));
    const strandLengthInput = Object.fromEntries(STRANDS.map(
      (strand) => {
        const input = ui.intInput(`${STRAND_NAME[strand]} length`, DEFAULT_SEQUENCE_LENGTH, () => updateUiForNewSequenceLength());
        input.setTooltip(`Length of ${STRAND_NAME[strand].toLowerCase()}, including overhangs`);
        return [strand, input];
      }));
    const strandVar = Object.fromEntries(STRANDS.map((strand) => [strand, '']));
    // const strandColumnInputDiv = Object.fromEntries(STRANDS.map(
    //   (strand) => [strand, ui.div([])]
    // ));
    const inputExample = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.textInput(
        ``, generateExample(strandLengthInput[strand].value!, sequenceBase.value!))
      ]));

    const strandColumnInput = Object.fromEntries(STRANDS.map((strand) => {
      const input = ui.choiceInput(`${STRAND_NAME[strand]} column`, '', [], (colName: string) => {
        validateStrandColumn(colName, strand);
        strandVar[strand] = colName;
      });
      // strandColumnInputDiv[strand].append(input.root);
      return [strand, input];
    }));

    const firstPto = Object.fromEntries(STRANDS.map((strand) => {
      const input = ui.boolInput(`First ${strand} PTO`, fullyPto.value!, () => updateSvgScheme());
      input.setTooltip(`ps linkage before first nucleotide of ${STRAND_NAME[strand].toLowerCase()}`);

      input.captionLabel.classList.add('ui-label-right');
      input.captionLabel.style.textAlign = 'left';
      input.captionLabel.style.maxWidth = '100px';
      input.captionLabel.style.minWidth = '40px';
      input.captionLabel.style.width = 'auto';

      return [strand, input];
    }));

    const terminalModification = Object.fromEntries(STRANDS.map((strand) => {
      const inputs = Object.fromEntries(TERMINAL_KEYS.map((key) => {
        const input = ui.stringInput(`${strand} ${TERMINAL[key]}\' Modification`, '', () => {
          updateSvgScheme();
          updateOutputExamples();
        });
        input.setTooltip(`Additional ${strand} ${TERMINAL[key]}\' Modification`);
        return [key, input];
      }));
      return [strand, inputs];
    }));

    const outputExample = Object.fromEntries(STRANDS.map((strand) => {
      const input = ui.textInput('', translateSequence(
        inputExample[strand].value, baseInputsObject[strand], ptoLinkages[strand], terminalModification[strand][THREE_PRIME],terminalModification[strand][FIVE_PRIME], firstPto[strand].value!
      ));
      input.input.style.minWidth = 'none';
      input.input.style.flexGrow = '1';
      $(input.root.lastChild).css('height', 'auto');
      return [strand, input];
    }));

    const modificationSection = Object.fromEntries(STRANDS.map((strand) => {
      const panel = ui.block([
        ui.h1(`${STRAND_NAME[strand]}`),
        ui.divH([
          ui.div([ui.divText('#')], {style: {width: '20px'}})!,
          ui.block75([ui.divText('Modification')])!,
          ui.div([ui.divText('PTO')])!,
        ]),
        modificationItems[strand],
      ], {style: {paddingTop: '12px'}});
      return [strand, panel];
    }));

    STRANDS.forEach((s) => {

      inputExample[s].input.style.resize = 'none';
      outputExample[s].input.style.resize = 'none';
      inputExample[s].input.style.minWidth = 'none';
      inputExample[s].input.style.flexGrow = '1';
      outputExample[s].input.style.minWidth = 'none';
      outputExample[s].input.style.flexGrow = '1';
      let options = ui.div([
        ui.button(ui.iconFA('copy', () => {}), () => {
          navigator.clipboard.writeText(outputExample[s].value).then(() =>
            grok.shell.info('Sequence was copied to clipboard'));
        }),
      ], 'ui-input-options');
      options.style.height = 'inherit';
      outputExample[s].root.append(
        options
      );
    })

    // const inputIdColumnDiv = ui.div([]);
    const svgDiv = ui.div([]);
    const asExampleDiv = ui.div([], 'ui-form ui-form-wide');
    const loadPatternDiv = ui.div([]);
    const asModificationDiv = ui.form([]);
    const isEnumerateModificationsDiv = ui.divH([
      ui.boolInput(defaultBase, true, (v: boolean) => {
        if (v) {
          if (!enumerateModifications.includes(defaultBase))
            enumerateModifications.push(defaultBase);
        } else {
          const index = enumerateModifications.indexOf(defaultBase, 0);
          if (index > -1)
            enumerateModifications.splice(index, 1);
        }
        updateSvgScheme();
        updateOutputExamples();
      }).root,
    ]);

    const asLengthDiv = ui.div([strandLengthInput[AS].root]);

    function getTableInput(tableList: DG.DataFrame[]): DG.InputBase<DG.DataFrame | null> {
      const tableInput = ui.tableInput('Tables', tableList[0], tableList, () => {
        const table = tableInput.value;
        if (table === null) {
          console.warn('Table is null');
          return;
        }
        const tableName = table!.name;
        if (!grok.shell.tableNames.includes(tableName)) {
          const view = grok.shell.v;
          grok.shell.addTableView(table!);
          grok.shell.v = view;
        }
        // console.log(`table:`, table);
        // console.log(`cols:`, table.columns);
        // console.log(`names:`, table.columns.names());
        const columnNames = table.columns.names();

        STRANDS.forEach((strand) => {
          const defaultColumn = columnNames[0];
          validateStrandColumn(defaultColumn, strand);
          strandVar[strand] = defaultColumn;
          const input = ui.choiceInput(`${STRAND_NAME[strand]} column`, defaultColumn, columnNames, (colName: string) => {
            validateStrandColumn(colName, strand);
            strandVar[strand] = colName;
            console.log(`clicked ${strand} var:`, strandVar[strand]);
          });
          $(strandColumnInput[strand].root).replaceWith(input.root);
          // $(inputIdColumnDiv).replaceWith(ui.div([inputStrandColumn[strand]]));
          // $(inputStrandColumnDiv[strand]).replaceWith(ui.divText('ababa'));
          // $(inputStrandColumnDiv[strand])
          //   .replaceWith(ui.div([inputStrandColumn[strand]]));
          // strandColumnInputDiv[strand].innerHTML = ui.divText('ababa').innerHTML;
          // // inputStrandColumnDiv[strand].append(inputStrandColumn[strand].root);
          // console.log(`${strand} input div`,strandColumnInputDiv[strand].innerHTML);
        })

        idVar = columnNames[0];
        // todo: unify with inputStrandColumn
        const idInput = ui.choiceInput('ID column', columnNames[0], columnNames, (colName: string) => {
          validateIdsColumn(colName);
          idVar = colName;
        });
        $(inputIdColumn.root).replaceWith(idInput.root);
      });
      return tableInput;
    }

    // const tableList = grok.shell.tables;
    const tableInput = getTableInput([]);

    // function updateInputs(): void {
    //   grok.shell.info('Event caught');
    //   const tableList = grok.shell.tables;
    //   console.log(`newTables:`, tableList);
    //   $(tableInput.root).replaceWith(getTableInput(tableList).root);
    // }

    // grok.events.onTableAdded.subscribe(() => updateInputs());
    // grok.events.onTableRemoved.subscribe(() => updateInputs());


    // todo: unify with strandVar
    let idVar = '';
    const inputIdColumn = ui.choiceInput('ID column', '', [], (colName: string) => {
      validateIdsColumn(colName);
      idVar = colName;
    });
    // inputIdColumnDiv.append(inputIdColumn.root);

    updatePatternsList();

    const createAsStrand = ui.boolInput('Anti sense strand', true, (v: boolean) => {
      modificationSection[AS].hidden = !v;
      // strandColumnInputDiv[AS].hidden = !v;
      strandColumnInput[AS].root.hidden = !v;
      asLengthDiv.hidden = !v;
      asModificationDiv.hidden = !v;
      asExampleDiv.hidden = !v;
      firstPto[AS].root.hidden = !v;
      updateSvgScheme();
    });
    createAsStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    const saveAs = ui.textInput('Save as', 'Pattern name', () => updateSvgScheme());
    saveAs.setTooltip('Name Of New Pattern');


    TERMINAL_KEYS.forEach((terminal) => {
      asModificationDiv.append(terminalModification[AS][terminal].root);
    })

    const comment = ui.textInput('Comment', '', () => updateSvgScheme());

    const savePatternButton = ui.bigButton('Save', () => {
      if (saveAs.value !== '')
        savePattern().then(() => grok.shell.info('Pattern saved'));
      else {
        const name = ui.stringInput('Enter name', '');
        ui.dialog('Pattern Name')
          .add(name.root)
          .onOK(() => {
            saveAs.value = name.value;
            savePattern().then(() => grok.shell.info('Pattern saved'));
          })
          .show();
      }
    });
    saveAs.addOptions(savePatternButton);

    const convertSequenceButton = ui.bigButton('Convert', () => {
      const condition = [true, createAsStrand.value];
      console.log(`strand vars:`, Object.values(strandVar));
      if (STRANDS.some((s, i) => condition[i] && strandVar[s] === ''))
        grok.shell.info('Please select table and columns on which to apply pattern');
      else if (STRANDS.some((s) => strandLengthInput[s].value !== inputExample[s].value.length)) {
        const dialog = ui.dialog('Length Mismatch');
        $(dialog.getButton('OK')).hide();
        dialog
          .add(ui.divText('Length of sequences in columns doesn\'t match entered length. Update length value?'))
          .addButton('YES', () => {
            STRANDS.forEach((s) => {
              strandLengthInput[s].value = tableInput.value!.getCol(strandColumnInput[s].value!).getString(0).length;
            })
            dialog.close();
          })
          .show();
      } else {
        if (idVar !== '')
          addColumnWithIds(tableInput.value!.name, idVar, getShortName(saveAs.value));
        const condition = [true, createAsStrand.value];
        STRANDS.forEach((strand, i) => {
          if (condition[i])
            addColumnWithTranslatedSequences(
              tableInput.value!.name, strandVar[strand], baseInputsObject[strand], ptoLinkages[strand],
              terminalModification[strand][FIVE_PRIME], terminalModification[strand][THREE_PRIME], firstPto[strand].value!);
        })
        grok.shell.v = grok.shell.getTableView(tableInput.value!.name);
        grok.shell.info(((createAsStrand.value) ? 'Columns were' : 'Column was') +
          ' added to table \'' + tableInput.value!.name + '\'');
        updateOutputExamples();
      }
    });

    asExampleDiv.append(inputExample[AS].root);
    asExampleDiv.append(outputExample[AS].root);

    updateUiForNewSequenceLength();

    const exampleSection = ui.div([
      ui.h1('Conversion preview'),
      inputExample[SS].root,
      outputExample[SS].root,
      asExampleDiv,
    ], 'ui-form ui-form-wide');

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

    return ui.splitH([
      ui.box(
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
          ui.buttonsInput([
            convertSequenceButton,
          ]),
        ], 'ui-form')
        , {style:{maxWidth:'450px'}}),
      ui.panel([
        svgDiv,
        isEnumerateModificationsDiv,
        ui.divH([
          downloadButton,
          editPattern
        ], {style:{gap:'12px', marginTop:'12px'}}),
        ui.divH([
          ui.divV([
            ui.h1('Sense strand'),
            inputExample[SS].root,
            outputExample[SS].root,
          ], 'ui-block'),
          ui.divV([
            ui.h1('Anti sense'),
            inputExample[AS],
            outputExample[AS]
          ], 'ui-block'),
        ], {style:{gap:'24px', marginTop:'24px'}}),
        ui.h1('Additional modifications'),
        ui.form([
          terminalModification[SS][FIVE_PRIME],
          terminalModification[SS][THREE_PRIME],
        ]),
        asModificationDiv,
      ], {style: {overflowX: 'scroll', padding:'12px 24px'}})
    ], {}, true)
  }
}
