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

type BooleanInput = DG.InputBase<boolean | null>;
type StringInput = DG.InputBase<string | null>;

export class PatternLayoutHandler {
  get htmlDivElement() {
    function updateModification(strand: string) {
      modificationItems[strand].innerHTML = '';
      ptoLinkages[strand] = ptoLinkages[strand].concat(Array(maxStrandLength[strand] - baseInputsObject[strand].length).fill(fullyPto));
      baseInputsObject[strand] = baseInputsObject[strand].concat(Array(maxStrandLength[strand] - baseInputsObject[strand].length).fill(sequenceBase));
      let nucleotideCounter = 0;
      for (let i = 0; i < strandLengthInput[strand].value!; i++) {
        ptoLinkages[strand][i] = ui.boolInput('', ptoLinkages[strand][i].value!, () => {
          updateSvgScheme();
          updateOutputExamples();
        });
        baseInputsObject[strand][i] = ui.choiceInput('', baseInputsObject[strand][i].value, baseChoices, (v: string) => {
          if (!enumerateModifications.includes(v)) {
            enumerateModifications.push(v);
            isEnumerateModificationsDiv.append(
              ui.divText('', {style: {width: '25px'}}),
              ui.boolInput(v, true, (boolV: boolean) => {
                if (boolV) {
                  if (!enumerateModifications.includes(v))
                    enumerateModifications.push(v);
                } else {
                  const index = enumerateModifications.indexOf(v, 0);
                  if (index > -1)
                    enumerateModifications.splice(index, 1);
                }
                updateSvgScheme();
              }).root,
            );
          }
          updateModification(AS);
          updateSvgScheme();
          updateOutputExamples();
        });
        $(baseInputsObject[strand][i].root).addClass('st-pattern-choice-input');
        if (!isOverhang(baseInputsObject[strand][i].value!))
          nucleotideCounter++;

        modificationItems[strand].append(
          ui.divH([
            ui.div([ui.label(isOverhang(baseInputsObject[strand][i].value!) ? '' : String(nucleotideCounter))],
              {style: {width: '20px'}})!,
            ui.block75([baseInputsObject[strand][i].root])!,
            ui.div([ptoLinkages[strand][i]])!,
          ], {style: {alignItems: 'center'}}),
        );
      }
    }

    function updateUiForNewSequenceLength() {
      if (Object.values(strandLengthInput).every((input) => input.value! < MAX_SEQUENCE_LENGTH)) {
        STRANDS.forEach((strand) => {
          if (strandLengthInput[strand].value! > maxStrandLength[strand])
            maxStrandLength[strand] = strandLengthInput[strand].value!;
          updateModification(strand);
        })

        updateSvgScheme();
        updateInputExamples();
        updateOutputExamples();
      } else {
        ui.dialog('Out of range')
          .add(ui.divText('Sequence length should be less than ' +
            MAX_SEQUENCE_LENGTH.toString() + ' due to UI constrains.'))
          .onOK(()=> {Object.values(strandLengthInput).every((input)=> input.value = 34)})
          .onCancel(()=> {Object.values(strandLengthInput).every((input)=> input.value = 34)})
          .showModal(false);
      }
    }

    // todo: unify with updateBases
    function updatePto(newPtoValue: boolean): void {
      STRANDS.forEach((strand) => {
        for (let i = 0; i < ptoLinkages[strand].length; i++)
          ptoLinkages[strand][i].value = newPtoValue;
      })
      updateSvgScheme();
    }

    function updateBases(newBasisValue: string): void {
      STRANDS.forEach((strand) => {
        for (let i = 0; i < baseInputsObject[strand].length; i++)
          baseInputsObject[strand][i].value = newBasisValue;
      })
      updateSvgScheme();
    }

    function updateInputExamples() {
      STRANDS.forEach((s) => {
        if (inputStrandColumn[s].value === '')
          inputExample[s].value = generateExample(strandLengthInput[s].value!, sequenceBase.value!);
      });
    }

    function updateOutputExamples() {
      const conditions = [true, createAsStrand.value];
      STRANDS.forEach((strand, i) => {
        if (conditions[i]) {
          outputExample[strand].value = translateSequence(inputExample[strand].value, baseInputsObject[strand], ptoLinkages[strand], terminalModification[strand][FIVE_PRIME], terminalModification[strand][THREE_PRIME], firstPto[strand].value!);
        }
      })
    }

    function updateSvgScheme() {
      svgDiv.innerHTML = '';
      svgDiv.append(
        ui.span([

          // todo: refactor the funciton, reduce # of args
          drawAxolabsPattern(
            getShortName(saveAs.value),
            createAsStrand.value!,

            baseInputsObject[SS].slice(0, strandLengthInput[SS].value!).map((e) => e.value!),
            baseInputsObject[AS].slice(0, strandLengthInput[AS].value!).map((e) => e.value!),

            [firstPto[SS].value!].concat(ptoLinkages[SS].slice(0, strandLengthInput[SS].value!).map((e) => e.value!)),
            [firstPto[AS].value!].concat(ptoLinkages[AS].slice(0, strandLengthInput[AS].value!).map((e) => e.value!)),

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
    function detectDefaultBasis(array: string[]) {
      const modeMap: {[index: string]: number} = {};
      let maxEl = array[0];
      let maxCount = 1;
      for (let i = 0; i < array.length; i++) {
        const el = array[i];
        if (modeMap[el] === null)
          modeMap[el] = 1;
        else
          modeMap[el]++;
        if (modeMap[el] > maxCount) {
          maxEl = el;
          maxCount = modeMap[el];
        }
      }
      return maxEl;
    }

    async function parsePatternAndUpdateUi(newName: string) {
      const pi = DG.TaskBarProgressIndicator.create('Loading pattern...');
      await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false).then((entities) => {
        const obj = JSON.parse(entities[newName]);
        sequenceBase.value = detectDefaultBasis(obj[FIELD.AS_BASES].concat(obj[FIELD.SS_BASES]));
        createAsStrand.value = (obj[FIELD.AS_BASES].length > 0);
        saveAs.value = newName;

        let fields = [FIELD.SS_BASES, FIELD.AS_BASES];
        STRANDS.forEach((strand, i) => {
          baseInputsObject[strand] = [];
          const field = fields[i];
          for (let j = 0; j < obj[field].length; j++)
            baseInputsObject[strand].push(ui.choiceInput('', obj[field][j], baseChoices));
        })

        fields = [FIELD.SS_PTO, FIELD.AS_PTO];
        STRANDS.forEach((s, i) => {
          const field = fields[i];
          firstPto[s].value = obj[field][0];
          ptoLinkages[s] = [];
          for (let j = 1; j < obj[field].length; j++)
            ptoLinkages[s].push(ui.boolInput('', obj[field][j]));
        });

        fields = [FIELD.SS_BASES, FIELD.AS_BASES];
        STRANDS.forEach((strand, i) => {
          strandLengthInput[strand].value = obj[fields[i]].length;
        })

        const field = [[FIELD.SS_3, FIELD.SS_5], [FIELD.AS_3, FIELD.AS_5]];
        STRANDS.forEach((strand, i) => {
          TERMINAL_KEYS.forEach((terminal, j) => {
            terminalModification[strand][terminal].value = obj[field[i][j]];
          })
        })
        comment.value = obj[FIELD.COMMENT];
      });
      pi.close();
    }

    function allColumnValuesOfEqualLength(colName: string): boolean {
      const col = tableInput.value!.getCol(colName);
      let allLengthsAreTheSame = true;
      for (let i = 1; i < col.length; i++) {
        if (col.get(i - 1).length !== col.get(i).length && col.get(i).length !== 0) {
          allLengthsAreTheSame = false;
          break;
        }
      }
      if (!allLengthsAreTheSame) {
        const dialog = ui.dialog('Sequences lengths mismatch');
        $(dialog.getButton('OK')).hide();
        dialog
          .add(ui.divText('The sequence length should match the number of Raw sequences in the input file'))
          .add(ui.divText('\'ADD COLUMN\' to see sequences lengths'))
          .addButton('ADD COLUMN', () => {
            tableInput.value!.columns.addNewInt('Sequences lengths in ' + colName).init((j: number) => col.get(j).length);
            grok.shell.info('Column with lengths added to \'' + tableInput.value!.name + '\'');
            dialog.close();
            grok.shell.v = grok.shell.getTableView(tableInput.value!.name);
          })
          .show();
      }
      if (col.get(0) !== strandLengthInput[SS].value) {
        const d = ui.dialog('Length was updated by value to from imported file');
        d.add(ui.divText('Latest modifications may not take effect during translation'))
          .onOK(() => grok.shell.info('Lengths changed')).show();
      }
      return allLengthsAreTheSame;
    }

    async function getCurrentUserName(): Promise<string> {
      return await grok.dapi.users.current().then((user) => {
        return ' (created by ' + user.friendlyName + ')';
      });
    }

    async function postPatternToUserStorage() {
      const currUserName = await getCurrentUserName();
      saveAs.value = (saveAs.stringValue.includes('(created by ')) ?
        getShortName(saveAs.value) + currUserName :
        saveAs.stringValue + currUserName;
      return grok.dapi.userDataStorage.postValue(
        USER_STORAGE_KEY,
        saveAs.value,
        JSON.stringify({
          [FIELD.SS_BASES]: baseInputsObject[SS].slice(0, strandLengthInput[SS].value!).map((e) => e.value),
          [FIELD.AS_BASES]: baseInputsObject[AS].slice(0, strandLengthInput[AS].value!).map((e) => e.value),
          [FIELD.SS_PTO]: [firstPto[SS].value].concat(ptoLinkages[SS].slice(0, strandLengthInput[SS].value!).map((e) => e.value)),
          [FIELD.AS_PTO]: [firstPto[AS].value].concat(ptoLinkages[AS].slice(0, strandLengthInput[AS].value!).map((e) => e.value)),
          [FIELD.SS_3]: terminalModification[SS][THREE_PRIME].value,
          [FIELD.SS_5]:terminalModification[SS][FIVE_PRIME].value,
          [FIELD.AS_3]: terminalModification[AS][THREE_PRIME].value,
          [FIELD.AS_5]: terminalModification[AS][FIVE_PRIME].value,
          [FIELD.COMMENT]: comment.value,
        }),
        false,
      ).then(() => grok.shell.info('Pattern \'' + saveAs.value + '\' was successfully uploaded!'));
    }

    async function updatePatternsList() {
      grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false).then(async (entities) => {
        const lstMy: string[] = [];
        const lstOthers: string[] = [];

        // TODO: display short name, but use long for querying userdataStorage
        for (const ent of Object.keys(entities)) {
          if (await isCurrentUserCreatedThisPattern(ent))
            lstOthers.push(ent);
          else
            lstMy.push(ent);//getShortName(ent));
        }

        let loadPattern = ui.choiceInput('Load pattern', '', lstMy, (v: string) => parsePatternAndUpdateUi(v));

        const currentUserName = (await grok.dapi.users.current()).friendlyName;
        const otherUsers = 'Other users';

        const patternListChoiceInput = ui.choiceInput('', currentUserName, [currentUserName, otherUsers], (v: string) => {
          const currentList = v === currentUserName ? lstMy : lstOthers;
          loadPattern = ui.choiceInput('Load pattern', '', currentList, (v: string) => parsePatternAndUpdateUi(v));

          loadPattern.root.append(patternListChoiceInput.input);
          loadPattern.root.append(loadPattern.input);
          // @ts-ignore
          loadPattern.input.style.maxWidth = '120px';
          loadPattern.input.style.marginLeft = '12px';
          loadPattern.setTooltip('Apply Existing Pattern');

          loadPatternDiv.innerHTML = '';
          loadPatternDiv.append(loadPattern.root);
          loadPattern.root.append(
            ui.div([
              ui.button(ui.iconFA('trash-alt', () => {}), async () => {
                if (loadPattern.value === null)
                  grok.shell.warning('Choose pattern to delete');
                else if (await isCurrentUserCreatedThisPattern(saveAs.value))
                  grok.shell.warning('Cannot delete pattern, created by other user');
                else {
                  await grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, loadPattern.value, false)
                    .then(() => grok.shell.info('Pattern \'' + loadPattern.value + '\' deleted'));
                }
                await updatePatternsList();
              }),
            ], 'ui-input-options'),
          );
        });
        patternListChoiceInput.input.style.maxWidth = '142px';
        loadPattern.root.append(patternListChoiceInput.input);
        loadPattern.root.append(loadPattern.input);
        // @ts-ignore
        loadPattern.input.style.maxWidth = '100px';
        loadPattern.setTooltip('Apply Existing Pattern');

        loadPatternDiv.innerHTML = '';
        loadPatternDiv.append(loadPattern.root);
        loadPattern.root.append(
          ui.div([
            ui.button(ui.iconFA('trash-alt', () => {}), async () => {
              if (loadPattern.value === null)
                grok.shell.warning('Choose pattern to delete');
              else if (await isCurrentUserCreatedThisPattern(saveAs.value))
                grok.shell.warning('Cannot delete pattern, created by other user');
              else {
                await grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, loadPattern.value, false)
                  .then(() => grok.shell.info('Pattern \'' + loadPattern.value + '\' deleted'));
              }
              await updatePatternsList();
            }),
          ], 'ui-input-options'),
        );
      });
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
      const allLengthsAreTheSame: boolean = allColumnValuesOfEqualLength(colName);
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
      updateBases(v);
      updateOutputExamples();
    });
    const fullyPto = ui.boolInput('Fully PTO', DEFAULT_PTO, (v: boolean) => {
      STRANDS.forEach((s) => { firstPto[s].value = v; })
      updatePto(v);
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
    // todo: rename to strandColumnInputDiv
    const inputStrandColumnDiv = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.div([])]
    ));
    const inputExample = Object.fromEntries(STRANDS.map(
      (strand) => [strand, ui.textInput(
        ``, generateExample(strandLengthInput[strand].value!, sequenceBase.value!))
      ]));

    // todo: rename to strandColumnInput
    const inputStrandColumn = Object.fromEntries(STRANDS.map((strand) => {
      const input: StringInput = ui.choiceInput(`${STRAND_NAME[strand]} column`, '', [], (colName: string) => {
        validateStrandColumn(colName, strand);
        strandVar[strand] = colName;
      });
      inputStrandColumnDiv[strand].append(input.root);
      //input.addPostfix(``);
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
      //inputExample[s].input.style.minWidth = EXAMPLE_MIN_WIDTH;
      outputExample[s].input.style.resize = 'none';
      //outputExample[s].input.styl
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

    const inputIdColumnDiv = ui.div([]);
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

    function getTableInput(tableList: DG.DataFrame[]): DG.InputBase {
      console.log(`tableList:`, tableList);
      const tableInput = ui.tableInput('Tables', tableList[0], tableList, (table: DG.DataFrame) => {
        console.log(`table:`, table);
        console.log(`cols:`, table.columns);
        console.log(`names:`, table.columns.names());

        STRANDS.forEach((strand) => {
          inputStrandColumn[strand] = ui.choiceInput(`${strand} column`, '', table.columns.names(), (colName: string) => {
            validateStrandColumn(colName, strand);
            strandVar[strand] = colName;
          });
          // $(inputIdColumnDiv).replaceWith(ui.div([inputStrandColumn[strand]]));
          inputStrandColumnDiv[strand].innerHTML = '';
          inputStrandColumnDiv[strand].append(inputStrandColumn[strand].root);
        })

        // todo: unify with inputStrandColumn
        console.log(`columns names:`, table);
        const inputIdColumn = ui.choiceInput('ID column', '', table.columns.names(), (colName: string) => {
          validateIdsColumn(colName);
          idVar = colName;
        });
        inputIdColumnDiv.innerHTML = '';
        inputIdColumnDiv.append(inputIdColumn.root);
      });
      return tableInput;
    }

    const tableList = grok.shell.tables;
    const tableInput = getTableInput(tableList);

    function updateInputs(): void {
      grok.shell.info('Event caught');
      const tableList = grok.shell.tables;
      console.log(`newTables:`, tableList);
      $(tableInput.root).replaceWith(getTableInput(tableList).root);
    }

    grok.events.onTableAdded.subscribe(() => updateInputs());
    grok.events.onTableRemoved.subscribe(() => updateInputs());


    // todo: unify with strandVar
    let idVar = '';
    const inputIdColumn = ui.choiceInput('ID column', '', [], (colName: string) => {
      validateIdsColumn(colName);
      idVar = colName;
    });
    inputIdColumnDiv.append(inputIdColumn.root);

    updatePatternsList();

    const createAsStrand = ui.boolInput('Anti sense strand', true, (v: boolean) => {
      modificationSection[AS].hidden = !v;
      inputStrandColumnDiv[AS].hidden = !v;
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
      if (STRANDS.some((s, i) => condition[i] && strandVar[s] === ''))
        grok.shell.info('Please select table and columns on which to apply pattern');
      else if (STRANDS.some((s) => strandLengthInput[s].value !== inputExample[s].value.length)) {
        const dialog = ui.dialog('Length Mismatch');
        $(dialog.getButton('OK')).hide();
        dialog
          .add(ui.divText('Length of sequences in columns doesn\'t match entered length. Update length value?'))
          .addButton('YES', () => {
            STRANDS.forEach((s) => {
              strandLengthInput[s].value = tableInput.value!.getCol(inputStrandColumn[s].value!).getString(0).length;
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
      inputStrandColumnDiv[SS],
      inputStrandColumnDiv[AS],
      inputIdColumnDiv,
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
          inputStrandColumn[SS],
          inputStrandColumn[AS],
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
