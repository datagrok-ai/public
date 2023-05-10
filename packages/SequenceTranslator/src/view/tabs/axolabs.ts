/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {JsonLoader, AxolabsStyle} from '../../model/data-loading-utils/json-loader';
import {defaultPto, defaultSequenceLength, maximalValidSequenceLength, userStorageKey, exampleMinWidth} from '../../model/axolabs/const';
import {isOverhang} from '../../model/axolabs/helpers';
import {generateExample, translateSequence, getShortName, isCurrentUserCreatedThisPattern, findDuplicates, addColumnWithIds, addColumnWithTranslatedSequences} from '../../model/axolabs/axolabs-tab';
import {drawAxolabsPattern} from '../../model/axolabs/draw-svg';
// todo: remove ts-ignore
//@ts-ignore
import * as svg from 'save-svg-as-png';

const enum FIELD {
  SS_BASES = 'ssBases',
  AS_BASES = 'asBases',
  SS_PTO = 'ssPtoLinkages',
  AS_PTO = 'asPtoLinkages',
  SS_3 = 'ssThreeModification',
  SS_5 = 'ssFiveModification',
  AS_3 = 'asThreeModification',
  AS_5 = 'asFiveModification',
  COMMENT = 'comment',
};

type BooleanInput = DG.InputBase<boolean | null>;     
type StringInput = DG.InputBase<string | null>;   
type NumberInput = DG.InputBase<number | null>;               

const IDX = {
  SS: 0,
  AS: 1,
  THREE_PRIME: 0,
  FIVE_PRIME: 1,
};

const strands = ['SS', 'AS'] as const;
const strandLongNames = ['Sense Strand', 'Antisense Strand'] as const;
const terminals = [3, 5] as const;

export class AxolabsTabUI {
  constructor() {
    this.axolabsStyle = JsonLoader.getInstance().getAxolabsStyleDictionary();
  }
  private axolabsStyle: AxolabsStyle;

  get htmlDivElement() {
    const baseChoices: string[] = Object.keys(this.axolabsStyle);
    const defaultBase: string = baseChoices[0];
    const enumerateModifications = [defaultBase];
    const sequenceBase = ui.choiceInput('Sequence Basis', defaultBase, baseChoices, (v: string) => {
      updateBases(v);
      updateOutputExamples();
    });
    const fullyPto = ui.boolInput('Fully PTO', defaultPto, (v: boolean) => {
      firstStrandPto[IDX.SS].value = v;
      firstStrandPto[IDX.AS].value = v;
      updatePto(v);
      updateOutputExamples();
    });



    const maximalStrandLength = strands.map(() => defaultSequenceLength);
    // todo: remove vague legacy 'items' from name
    const strandModificationItems = strands.map(() => ui.div([]));
    const strandPtoLinkages = strands.map(() =>
      Array<BooleanInput>(defaultSequenceLength).fill(ui.boolInput('', defaultPto))
    );
    const strandBases = strands.map(() => 
      Array<StringInput>(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices))
    );
    const strandLengthInput = strands.map((strand, i) => {
      const input = ui.intInput(`${strand} Length`, defaultSequenceLength, () => updateUiForNewSequenceLength());
      input.setTooltip(`Length of ${strandLongNames[i].toLowerCase()}, including overhangs`);
      return input;
    })
    const strandVar = strands.map(() => '');
    // todo: rename to strandColumnInputDiv
    const inputStrandColumnDiv = strands.map(() => ui.div([]));
    const strandInputExample = strands.map((_, i) => {
      return ui.textInput(`${strandLongNames[i]}`, generateExample(strandLengthInput[i].value!, sequenceBase.value!));
    })
   
    // todo: rename to strandColumnInput
    const inputStrandColumn = strands.map((_, i) => {
      const input: StringInput = ui.choiceInput(`${strandLongNames[i]} Column`, '', [], (colName: string) => {
        validateStrandColumn(colName, IDX.SS);
        strandVar[i] = colName;
      });
      inputStrandColumnDiv[i].append(input.root);
      return input;
    })

    const firstStrandPto = strands.map((strand, i) => {
      const input = ui.boolInput(`First ${strand} PTO`, fullyPto.value!, () => updateSvgScheme());
      input.setTooltip(`ps linkage before first nucleotide of ${strandLongNames[i].toLowerCase()}`);
      return input;
    });

    const strandTerminalModification = strands.map((strand) => {
      return terminals.map((terminal) => {
        const input = ui.stringInput(`${strand} ${terminal}\' Modification`, '', () => {
          updateSvgScheme();
          updateOutputExamples();
        });

        input.setTooltip(`Additional ${strand} ${terminal}\' Modification`);
        return input;
      })
    });

    const strandOutputExample = strands.map((_, i) =>
      ui.textInput(' ', translateSequence(
        strandInputExample[i].value, strandBases[i], strandPtoLinkages[i], strandTerminalModification[i][IDX.THREE_PRIME],strandTerminalModification[i][IDX.FIVE_PRIME], firstStrandPto[i].value!
      ))
    );

    strandInputExample[IDX.SS].input.style.resize = 'none';
    strandInputExample[IDX.AS].input.style.resize = 'none';

    strandInputExample[IDX.SS].input.style.minWidth = exampleMinWidth;
    strandInputExample[IDX.AS].input.style.minWidth = exampleMinWidth;

    strandOutputExample[IDX.SS].input.style.resize = 'none';
    strandOutputExample[IDX.AS].input.style.resize = 'none';

    strandOutputExample[IDX.SS].input.style.minWidth = exampleMinWidth;
    strandOutputExample[IDX.AS].input.style.minWidth = exampleMinWidth;

    // @ts-ignore
    strandOutputExample[IDX.SS].input.disabled = 'true';
    // @ts-ignore
    strandOutputExample[IDX.AS].input.disabled = 'true';

    strandOutputExample[IDX.SS].root.append(
      ui.div([
        ui.button(ui.iconFA('copy', () => {}), () => {
          navigator.clipboard.writeText(strandOutputExample[IDX.SS].value).then(() =>
            grok.shell.info('Sequence was copied to clipboard'));
        }),
      ], 'ui-input-options'),
    );
    strandOutputExample[IDX.AS].root.append(
      ui.div([
        ui.button(ui.iconFA('copy', () => {}), () => {
          navigator.clipboard.writeText(strandOutputExample[IDX.AS].value).then(() =>
            grok.shell.info('Sequence was copied to clipboard'));
        }),
      ], 'ui-input-options'),
    );


    const firstAsPtoDiv = ui.div([]);
    firstAsPtoDiv.append(firstStrandPto[IDX.AS].root);


    function updateAsModification() {
      strandModificationItems[IDX.AS].innerHTML = '';
      strandPtoLinkages[IDX.AS] = strandPtoLinkages[IDX.AS].concat(Array(maximalStrandLength[IDX.AS] - strandBases[IDX.AS].length).fill(fullyPto));
      strandBases[IDX.AS] = strandBases[IDX.AS].concat(Array(maximalStrandLength[IDX.AS] - strandBases[IDX.AS].length).fill(sequenceBase));
      let nucleotideCounter = 0;
      for (let i = 0; i < strandLengthInput[IDX.AS].value!; i++) {
        strandPtoLinkages[IDX.AS][i] = ui.boolInput('', strandPtoLinkages[IDX.AS][i].value!, () => {
          updateSvgScheme();
          updateOutputExamples();
        });
        strandBases[IDX.AS][i] = ui.choiceInput('', strandBases[IDX.AS][i].value, baseChoices, (v: string) => {
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
          updateAsModification();
          updateSvgScheme();
          updateOutputExamples();
        });
        if (!isOverhang(strandBases[IDX.AS][i].value!))
          nucleotideCounter++;

        strandModificationItems[IDX.AS].append(
          ui.divH([
            ui.div([ui.label(isOverhang(strandBases[IDX.AS][i].value!) ? '' : String(nucleotideCounter))],
              {style: {width: '20px'}})!,
            ui.block75([strandBases[IDX.AS][i].root])!,
            ui.div([strandPtoLinkages[IDX.AS][i]])!,
          ], {style: {alignItems: 'center'}}),
        );
      }
    }

    function updateSsModification() {
      strandModificationItems[IDX.SS].innerHTML = '';
      strandPtoLinkages[IDX.SS] = strandPtoLinkages[IDX.SS].concat(Array(maximalStrandLength[IDX.SS] - strandBases[IDX.SS].length).fill(fullyPto));
      strandBases[IDX.SS] = strandBases[IDX.SS].concat(Array(maximalStrandLength[IDX.SS] - strandBases[IDX.SS].length).fill(sequenceBase));
      let nucleotideCounter = 0;
      for (let i = 0; i < strandLengthInput[IDX.SS].value!; i++) {
        strandPtoLinkages[IDX.SS][i] = ui.boolInput('', strandPtoLinkages[IDX.SS][i].value!, () => {
          updateSvgScheme();
          updateOutputExamples();
        });
        strandBases[IDX.SS][i] = ui.choiceInput('', strandBases[IDX.SS][i].value, baseChoices, (v: string) => {
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
          updateSsModification();
          updateSvgScheme();
          updateOutputExamples();
        });
        if (!isOverhang(strandBases[IDX.SS][i].value!))
          nucleotideCounter++;

        strandModificationItems[IDX.SS].append(
          ui.divH([
            ui.div([ui.label(isOverhang(strandBases[IDX.SS][i].value!) ? '' : String(nucleotideCounter))],
              {style: {width: '20px'}})!,
            ui.block75([strandBases[IDX.SS][i].root])!,
            ui.div([strandPtoLinkages[IDX.SS][i]])!,
          ], {style: {alignItems: 'center'}}),
        );
      }
    }

    function updateUiForNewSequenceLength() {
      if (strandLengthInput[IDX.SS].value! < maximalValidSequenceLength && strandLengthInput[IDX.AS].value! < maximalValidSequenceLength) {
        if (strandLengthInput[IDX.SS].value! > maximalStrandLength[IDX.SS])
          maximalStrandLength[IDX.SS] = strandLengthInput[IDX.SS].value!;
        if (strandLengthInput[IDX.AS].value! > maximalStrandLength[IDX.AS])
          maximalStrandLength[IDX.AS] = strandLengthInput[IDX.AS].value!;
        updateSsModification();
        updateAsModification();
        updateSvgScheme();
        updateInputExamples();
        updateOutputExamples();
      } else {
        ui.dialog('Sequence length is out of range')
          .add(ui.divText('Sequence length should be less than ' +
            maximalValidSequenceLength.toString() + ' due to UI constrains.'))
          .add(ui.divText('Please change sequence length in order to define new pattern.'))
          .show();
      }
    }

    function updatePto(newPtoValue: boolean): void {
      for (let i = 0; i < strandPtoLinkages[IDX.SS].length; i++)
      strandPtoLinkages[IDX.SS][i].value = newPtoValue;

      for (let i = 0; i < strandPtoLinkages[IDX.AS].length; i++)
      strandPtoLinkages[IDX.AS][i].value = newPtoValue;

      updateSvgScheme();
    }

    function updateBases(newBasisValue: string): void {
      for (let i = 0; i < strandBases[IDX.SS].length; i++)
      strandBases[IDX.SS][i].value = newBasisValue;

      for (let i = 0; i < strandBases[IDX.AS].length; i++)
      strandBases[IDX.AS][i].value = newBasisValue;

      updateSvgScheme();
    }

    function updateInputExamples() {
      if (inputStrandColumn[IDX.SS].value == '')
        strandInputExample[IDX.SS].value = generateExample(strandLengthInput[IDX.SS].value!, sequenceBase.value!);
      if (createAsStrand.value && inputStrandColumn[IDX.AS].value == '')
        strandInputExample[IDX.AS].value = generateExample(strandLengthInput[IDX.AS].value!, sequenceBase.value!);
    }

    function updateOutputExamples() {
      strandOutputExample[IDX.SS].value = translateSequence(
        strandInputExample[IDX.SS].value, strandBases[IDX.SS], strandPtoLinkages[IDX.SS],strandTerminalModification[IDX.SS][IDX.FIVE_PRIME], strandTerminalModification[IDX.SS][IDX.THREE_PRIME], firstStrandPto[IDX.SS].value!);
      if (createAsStrand.value) {
        strandOutputExample[IDX.AS].value = translateSequence(
          strandInputExample[IDX.AS].value, strandBases[IDX.AS], strandPtoLinkages[IDX.AS], strandTerminalModification[IDX.AS][IDX.FIVE_PRIME], strandTerminalModification[IDX.AS][IDX.THREE_PRIME], firstStrandPto[IDX.AS].value!);
      }
    }

    function updateSvgScheme() {
      svgDiv.innerHTML = '';
      svgDiv.append(
        ui.span([
          drawAxolabsPattern(getShortName(saveAs.value),
            createAsStrand.value!,
            strandBases[IDX.SS].slice(0, strandLengthInput[IDX.SS].value!).map((e) => e.value!),
            strandBases[IDX.AS].slice(0, strandLengthInput[IDX.AS].value!).map((e) => e.value!),
            [firstStrandPto[IDX.SS].value!].concat(strandPtoLinkages[IDX.SS].slice(0, strandLengthInput[IDX.SS].value!).map((e) => e.value!)),
            [firstStrandPto[IDX.AS].value!].concat(strandPtoLinkages[IDX.AS].slice(0, strandLengthInput[IDX.AS].value!).map((e) => e.value!)),
            strandTerminalModification[IDX.SS][IDX.THREE_PRIME].value,
           strandTerminalModification[IDX.SS][IDX.FIVE_PRIME].value,
            strandTerminalModification[IDX.AS][IDX.THREE_PRIME].value,
            strandTerminalModification[IDX.AS][IDX.FIVE_PRIME].value,
            comment.value,
            enumerateModifications,
          ),
        ]),
      );
    }

    function detectDefaultBasis(array: string[]) {
      const modeMap: {[index: string]: number} = {};
      let maxEl = array[0];
      let maxCount = 1;
      for (let i = 0; i < array.length; i++) {
        const el = array[i];
        if (modeMap[el] == null)
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
      await grok.dapi.userDataStorage.get(userStorageKey, false).then((entities) => {
        const obj = JSON.parse(entities[newName]);
        sequenceBase.value = detectDefaultBasis(obj[FIELD.AS_BASES].concat(obj[FIELD.SS_BASES]));
        createAsStrand.value = (obj[FIELD.AS_BASES].length > 0);
        saveAs.value = newName;

        strandBases[IDX.SS] = [];
        for (let i = 0; i < obj[FIELD.SS_BASES].length; i++)
          strandBases[IDX.SS].push(ui.choiceInput('', obj[FIELD.SS_BASES][i], baseChoices));

        strandBases[IDX.AS] = [];
        for (let i = 0; i < obj[FIELD.AS_BASES].length; i++)
          strandBases[IDX.AS].push(ui.choiceInput('', obj[FIELD.AS_BASES][i], baseChoices));

        firstStrandPto[IDX.SS].value = obj[FIELD.SS_PTO][0];
        strandPtoLinkages[IDX.SS] = [];
        for (let i = 1; i < obj[FIELD.SS_PTO].length; i++)
          strandPtoLinkages[IDX.SS].push(ui.boolInput('', obj[FIELD.SS_PTO][i]));

        firstStrandPto[IDX.AS].value = obj[FIELD.AS_PTO][0];
        strandPtoLinkages[IDX.AS] = [];
        for (let i = 1; i < obj[FIELD.AS_PTO].length; i++)
          strandPtoLinkages[IDX.AS].push(ui.boolInput('', obj[FIELD.AS_PTO][i]));

        strandLengthInput[IDX.SS].value = obj[FIELD.SS_BASES].length;
        strandLengthInput[IDX.AS].value = obj[FIELD.AS_BASES].length;

        strandTerminalModification[IDX.SS][IDX.THREE_PRIME].value = obj[FIELD.SS_3];
       strandTerminalModification[IDX.SS][IDX.FIVE_PRIME].value = obj[FIELD.SS_5];
        strandTerminalModification[IDX.AS][IDX.THREE_PRIME].value = obj[FIELD.AS_3];
        strandTerminalModification[IDX.AS][IDX.FIVE_PRIME].value = obj[FIELD.AS_5];
        comment.value = obj[FIELD.COMMENT];
      });
      pi.close();
    }

    function checkWhetherAllValuesInColumnHaveTheSameLength(colName: string): boolean {
      const col = tables.value!.getCol(colName);
      let allLengthsAreTheSame = true;
      for (let i = 1; i < col.length; i++) {
        if (col.get(i - 1).length != col.get(i).length && col.get(i).length != 0) {
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
            tables.value!.columns.addNewInt('Sequences lengths in ' + colName).init((j: number) => col.get(j).length);
            grok.shell.info('Column with lengths added to \'' + tables.value!.name + '\'');
            dialog.close();
            grok.shell.v = grok.shell.getTableView(tables.value!.name);
          })
          .show();
      }
      if (col.get(0) != strandLengthInput[IDX.SS].value) {
        const d = ui.dialog('Length was updated by value to from imported file');
        d.add(ui.divText('Latest modifications may not take effect during translation'))
          .onOK(() => grok.shell.info('Lengths changed')).show();
      }
      return allLengthsAreTheSame;
    }

    async function getCurrentUserName(): Promise<string> {
      return await grok.dapi.users.current().then((user) => {
        return ' (created by ' + user.firstName + ' ' + user.lastName + ')';
      });
    }

    async function postPatternToUserStorage() {
      const currUserName = await getCurrentUserName();
      saveAs.value = (saveAs.stringValue.includes('(created by ')) ?
        getShortName(saveAs.value) + currUserName :
        saveAs.stringValue + currUserName;
      return grok.dapi.userDataStorage.postValue(
        userStorageKey,
        saveAs.value,
        JSON.stringify({
          [FIELD.SS_BASES]: strandBases[IDX.SS].slice(0, strandLengthInput[IDX.SS].value!).map((e) => e.value),
          [FIELD.AS_BASES]: strandBases[IDX.AS].slice(0, strandLengthInput[IDX.AS].value!).map((e) => e.value),
          [FIELD.SS_PTO]: [firstStrandPto[IDX.SS].value].concat(strandPtoLinkages[IDX.SS].slice(0, strandLengthInput[IDX.SS].value!).map((e) => e.value)),
          [FIELD.AS_PTO]: [firstStrandPto[IDX.AS].value].concat(strandPtoLinkages[IDX.AS].slice(0, strandLengthInput[IDX.AS].value!).map((e) => e.value)),
          [FIELD.SS_3]: strandTerminalModification[IDX.SS][IDX.THREE_PRIME].value,
          [FIELD.SS_5]:strandTerminalModification[IDX.SS][IDX.FIVE_PRIME].value,
          [FIELD.AS_3]: strandTerminalModification[IDX.AS][IDX.THREE_PRIME].value,
          [FIELD.AS_5]: strandTerminalModification[IDX.AS][IDX.FIVE_PRIME].value,
          [FIELD.COMMENT]: comment.value,
        }),
        false,
      ).then(() => grok.shell.info('Pattern \'' + saveAs.value + '\' was successfully uploaded!'));
    }

    async function updatePatternsList() {
      grok.dapi.userDataStorage.get(userStorageKey, false).then(async (entities) => {
        const lstMy: string[] = [];
        const lstOthers: string[] = [];

        // TODO: display short name, but use long for querying userdataStorage
        for (const ent of Object.keys(entities)) {
          if (await isCurrentUserCreatedThisPattern(ent))
            lstOthers.push(ent);
          else
            lstMy.push(ent);//getShortName(ent));
        }

        let loadPattern = ui.choiceInput('Load Pattern', '', lstMy, (v: string) => parsePatternAndUpdateUi(v));

        const myOrOthersPatternList = ui.choiceInput('', 'Mine', ['Mine', 'Others'], (v: string) => {
          const currentList = v == 'Mine' ? lstMy : lstOthers;
          loadPattern = ui.choiceInput('Load Pattern', '', currentList, (v: string) => parsePatternAndUpdateUi(v));

          loadPattern.root.append(myOrOthersPatternList.input);
          loadPattern.root.append(loadPattern.input);
          // @ts-ignore
          loadPattern.input.style.maxWidth = '100px';
          loadPattern.setTooltip('Apply Existing Pattern');

          loadPatternDiv.innerHTML = '';
          loadPatternDiv.append(loadPattern.root);
          loadPattern.root.append(
            ui.div([
              ui.button(ui.iconFA('trash-alt', () => {}), async () => {
                if (loadPattern.value == null)
                  grok.shell.warning('Choose pattern to delete');
                else if (await isCurrentUserCreatedThisPattern(saveAs.value))
                  grok.shell.warning('Cannot delete pattern, created by other user');
                else {
                  await grok.dapi.userDataStorage.remove(userStorageKey, loadPattern.value, false)
                    .then(() => grok.shell.info('Pattern \'' + loadPattern.value + '\' deleted'));
                }
                await updatePatternsList();
              }),
            ], 'ui-input-options'),
          );
        });
        loadPattern.root.append(myOrOthersPatternList.input);
        loadPattern.root.append(loadPattern.input);
        // @ts-ignore
        loadPattern.input.style.maxWidth = '100px';
        loadPattern.setTooltip('Apply Existing Pattern');

        loadPatternDiv.innerHTML = '';
        loadPatternDiv.append(loadPattern.root);
        loadPattern.root.append(
          ui.div([
            ui.button(ui.iconFA('trash-alt', () => {}), async () => {
              if (loadPattern.value == null)
                grok.shell.warning('Choose pattern to delete');
              else if (await isCurrentUserCreatedThisPattern(saveAs.value))
                grok.shell.warning('Cannot delete pattern, created by other user');
              else {
                await grok.dapi.userDataStorage.remove(userStorageKey, loadPattern.value, false)
                  .then(() => grok.shell.info('Pattern \'' + loadPattern.value + '\' deleted'));
              }
              await updatePatternsList();
            }),
          ], 'ui-input-options'),
        );
      });
    }

    async function savePattern() {
      await grok.dapi.userDataStorage.get(userStorageKey, false)
        .then((entities) => {
          if (Object.keys(entities).includes(saveAs.value)) {
            const dialog = ui.dialog('Pattern already exists');
            $(dialog.getButton('OK')).hide();
            dialog
              .add(ui.divText('Pattern name \'' + saveAs.value + '\' already exists.'))
              .add(ui.divText('Replace pattern?'))
              .addButton('YES', async () => {
                await grok.dapi.userDataStorage.remove(userStorageKey, saveAs.value, false)
                  .then(() => postPatternToUserStorage());
                dialog.close();
              })
              .show();
          } else
            postPatternToUserStorage();
        });
      await updatePatternsList();
    }

    const inputIdColumnDiv = ui.div([]);
    const svgDiv = ui.div([]);
    const asExampleDiv = ui.div([]);
    const appAxolabsDescription = ui.div([]);
    const loadPatternDiv = ui.div([]);
    const asModificationDiv = ui.div([]);
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

    const asLengthDiv = ui.div([strandLengthInput[IDX.AS].root]);
    
    function validateStrandColumn(colName: string, strandIdx: number): void {
      const allLengthsAreTheSame: boolean = checkWhetherAllValuesInColumnHaveTheSameLength(colName);
      const firstSequence = tables.value!.getCol(colName).get(0);
      if (allLengthsAreTheSame && firstSequence.length != strandLengthInput[strandIdx].value)
      strandLengthInput[strandIdx].value = tables.value!.getCol(colName).get(0).length;
      strandInputExample[strandIdx].value = firstSequence;
    }

    function validateIdsColumn(colName: string) {
      const col = tables.value!.getCol(colName);
      if (col.type != DG.TYPE.INT)
        grok.shell.error('Column should contain integers only');
      else if (col.categories.filter((e) => e != '').length < col.toList().filter((e) => e != '').length) {
        const duplicates = findDuplicates(col.getRawData());
        ui.dialog('Non-unique IDs')
          .add(ui.divText('Press \'OK\' to select rows with non-unique values'))
          .onOK(() => {
            const selection = tables.value!.selection;
            selection.init((i: number) => duplicates.indexOf(col.get(i)) > -1);
            grok.shell.v = grok.shell.getTableView(tables.value!.name);
            grok.shell.info('Rows are selected in table \'' + tables.value!.name + '\'');
          })
          .show();
      }
    }

    const tables = ui.tableInput('Tables', grok.shell.tables[0], grok.shell.tables, (t: DG.DataFrame) => {
      inputStrandColumn[IDX.SS] = ui.choiceInput('SS Column', '', t.columns.names(), (colName: string) => {
        validateStrandColumn(colName, IDX.SS);
        strandVar[IDX.SS] = colName;
      });
      inputStrandColumnDiv[IDX.SS].innerHTML = '';
      inputStrandColumnDiv[IDX.SS].append(inputStrandColumn[IDX.SS].root);
      inputStrandColumn[IDX.AS] = ui.choiceInput('AS Column', '', t.columns.names(), (colName: string) => {
        validateStrandColumn(colName, IDX.AS);
        strandVar[IDX.AS] = colName;
      });
      inputStrandColumnDiv[IDX.AS].innerHTML = '';
      inputStrandColumnDiv[IDX.AS].append(inputStrandColumn[IDX.AS].root);
      const inputIdColumn = ui.choiceInput('ID Column', '', t.columns.names(), (colName: string) => {
        validateIdsColumn(colName);
        idVar = colName;
      });
      inputIdColumnDiv.innerHTML = '';
      inputIdColumnDiv.append(inputIdColumn.root);
    });


    let idVar = '';
    const inputIdColumn = ui.choiceInput('ID Column', '', [], (colName: string) => {
      validateIdsColumn(colName);
      idVar = colName;
    });
    inputIdColumnDiv.append(inputIdColumn.root);

    updatePatternsList();

    const createAsStrand = ui.boolInput('Create AS Strand', true, (v: boolean) => {
      asModificationSection.hidden = (!v);
      inputStrandColumnDiv[IDX.AS].hidden = (!v);
      asLengthDiv.hidden = (!v);
      asModificationDiv.hidden = (!v);
      asExampleDiv.hidden = (!v);
      firstAsPtoDiv.hidden = (!v);
      updateSvgScheme();
    });
    createAsStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    const saveAs = ui.textInput('Save As', 'Pattern Name', () => updateSvgScheme());
    saveAs.setTooltip('Name Of New Pattern');


    asModificationDiv.append(strandTerminalModification[IDX.AS][IDX.THREE_PRIME].root);
    asModificationDiv.append(strandTerminalModification[IDX.AS][IDX.FIVE_PRIME].root);

    const comment = ui.textInput('Comment', '', () => updateSvgScheme());

    const savePatternButton = ui.button('Save', () => {
      if (saveAs.value != '')
        savePattern().then((r) => grok.shell.info('Pattern saved'));
      else {
        const name = ui.stringInput('Enter Name', '');
        ui.dialog('Pattern Name')
          .add(name.root)
          .onOK(() => {
            saveAs.value = name.value;
            savePattern().then((r) => grok.shell.info('Pattern saved'));
          })
          .show();
      }
    });

    const convertSequenceButton = ui.button('Convert Sequences', () => {
      if (strandVar[IDX.SS] == '' || (createAsStrand.value && strandVar[IDX.AS] == ''))
        grok.shell.info('Please select table and columns on which to apply pattern');
      else if (strandLengthInput[IDX.SS].value != strandInputExample[IDX.SS].value.length || strandLengthInput[IDX.AS].value != strandInputExample[IDX.AS].value.length) {
        const dialog = ui.dialog('Length Mismatch');
        $(dialog.getButton('OK')).hide();
        dialog
          .add(ui.divText('Length of sequences in columns doesn\'t match entered length. Update length value?'))
          .addButton('YES', () => {
            strandLengthInput[IDX.SS].value = tables.value!.getCol(inputStrandColumn[IDX.SS].value!).getString(0).length;
            strandLengthInput[IDX.AS].value = tables.value!.getCol(inputStrandColumn[IDX.AS].value!).getString(0).length;
            dialog.close();
          })
          .show();
      } else {
        if (idVar != '')
          addColumnWithIds(tables.value!.name, idVar, getShortName(saveAs.value));
        addColumnWithTranslatedSequences(
          tables.value!.name, strandVar[IDX.SS], strandBases[IDX.SS], strandPtoLinkages[IDX.SS],
         strandTerminalModification[IDX.SS][IDX.FIVE_PRIME], strandTerminalModification[IDX.SS][IDX.THREE_PRIME], firstStrandPto[IDX.SS].value!);
        if (createAsStrand.value) {
          addColumnWithTranslatedSequences(
            tables.value!.name, strandVar[IDX.AS], strandBases[IDX.AS], strandPtoLinkages[IDX.AS],
            strandTerminalModification[IDX.AS][IDX.FIVE_PRIME], strandTerminalModification[IDX.AS][IDX.THREE_PRIME], firstStrandPto[IDX.AS].value!);
        }
        grok.shell.v = grok.shell.getTableView(tables.value!.name);
        grok.shell.info(((createAsStrand.value) ? 'Columns were' : 'Column was') +
          ' added to table \'' + tables.value!.name + '\'');
        updateOutputExamples();
      }
    });

    asExampleDiv.append(strandInputExample[IDX.AS].root);
    asExampleDiv.append(strandOutputExample[IDX.AS].root);

    updateUiForNewSequenceLength();

    const exampleSection = ui.div([
      ui.h1('Example'),
      strandInputExample[IDX.SS].root,
      strandOutputExample[IDX.SS].root,
      asExampleDiv,
    ], 'ui-form');

    const inputsSection = ui.div([
      ui.h1('Inputs'),
      ui.divH([
        tables.root,
        inputStrandColumnDiv[IDX.SS],
      ]),
      ui.divH([
        inputStrandColumnDiv[IDX.AS],
        inputIdColumnDiv,
      ]),
      ui.buttonsInput([
        convertSequenceButton,
      ]),
    ], 'ui-form');

    const downloadButton = ui.button('Download', () => svg.saveSvgAsPng(document.getElementById('mySvg'), saveAs.value,
      {backgroundColor: 'white'}));

    const mainSection = ui.panel([
      ui.block([
        svgDiv,
      ], {style: {overflowX: 'scroll'}}),
      downloadButton,
      isEnumerateModificationsDiv,
      ui.div([
        ui.div([
          ui.divH([
            ui.h1('Pattern'),
            ui.div([
              ui.iconFA('question-circle', () => {
                appAxolabsDescription.innerHTML = '';
                appAxolabsDescription.append(info);
              }),
            ], {style: {padding: '2px'}}),
          ]),
          ui.divH([
            ui.div([
              strandLengthInput[IDX.SS].root,
              asLengthDiv,
              sequenceBase.root,
              comment.root,
              loadPatternDiv,
              saveAs.root,
              ui.buttonsInput([
                savePatternButton,
              ]),
            ], 'ui-form'),
            ui.div([
              createAsStrand.root,
              fullyPto.root,
              firstStrandPto[IDX.SS].root,
              firstAsPtoDiv,
              strandTerminalModification[IDX.SS][IDX.FIVE_PRIME].root,
              strandTerminalModification[IDX.SS][IDX.THREE_PRIME].root,
              asModificationDiv,
            ], 'ui-form'),
          ], 'ui-form'),
        ], 'ui-form'),
        inputsSection,
        exampleSection,
      ], {style: {flexWrap: 'wrap'}}),
    ]);

    const ssModificationSection = ui.panel([
      ui.h1('Sense Strand'),
      ui.divH([
        ui.div([ui.divText('#')], {style: {width: '20px'}})!,
        ui.block75([ui.divText('Modification')])!,
        ui.div([ui.divText('PTO')], {style: {paddingRight: '8px'}})!,
      ]),
      strandModificationItems[IDX.SS],
    ])!;

    const asModificationSection = ui.panel([
      ui.h1('Antisense Strand'),
      ui.divH([
        ui.div([ui.divText('#')], {style: {width: '20px'}})!,
        ui.block75([ui.divText('Modification')])!,
        ui.div([ui.divText('PTO')], {style: {paddingRight: '8px'}})!,
      ]),
      strandModificationItems[IDX.AS],
    ])!;

    const info = ui.info(
      [
        ui.divText('\n How to define new pattern:', {style: {'font-weight': 'bolder'}}),
        ui.divText('1. Choose table and columns with sense and antisense strands'),
        ui.divText('2. Choose lengths of both strands by editing checkboxes below'),
        ui.divText('3. Choose basis and PTO status for each nucleotide'),
        ui.divText('4. Set additional modifications for sequence edges'),
        ui.divText('5. Press \'Convert Sequences\' button'),
        ui.divText('This will add the result column(s) to the right of the table'),
      ], 'Create and apply Axolabs translation patterns.',
    );

    return ui.splitH([
      ui.div([
        appAxolabsDescription,
        mainSection!,
      ])!,
      ui.box(
        ui.divH([
          ssModificationSection,
          asModificationSection,
        ]), {style: {maxWidth: '360px'}},
      ),
    ]);
  }
}
