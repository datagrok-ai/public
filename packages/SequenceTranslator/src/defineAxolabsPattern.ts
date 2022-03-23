/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// @ts-ignore
import * as svg from 'save-svg-as-png';
import $ from 'cash-dom';

import {drawAxolabsPattern} from './drawAxolabsPattern';
import {axolabsMap} from './axolabsMap';

const baseChoices: string[] = Object.keys(axolabsMap);
const defaultBase: string = baseChoices[0];
const defaultPto: boolean = true;
const defaultSequenceLength: number = 23;
const maximalValidSequenceLength: number = 35;
const userStorageKey: string = 'SequenceTranslator';
const exampleMinWidth: string = '400px';

function generateExample(sequenceLength: number, sequenceBasis: string): string {
  const uniqueSymbols = axolabsMap[sequenceBasis].symbols.join('');
  return uniqueSymbols.repeat(Math.floor(sequenceLength / 4)) + uniqueSymbols.slice(0, sequenceLength % 4);
}

function findDuplicates(data: Int32Array | Float32Array | Float64Array | Uint32Array): number[] {
  return Array.from(new Set(data)).filter((value) => data.indexOf(value) !== data.lastIndexOf(value));
}

async function isCurrentUserCreatedThisPattern(patternName: string): Promise<boolean> {
  return await grok.dapi.users.current().then((user) => {
    const [firstName, lastName] = getUserName(patternName);
    return (user.firstName != firstName || user.lastName != lastName);
  });
}

function getShortName(patternName: string): string {
  let first = patternName.length + 1;
  for (let i = 0; i < patternName.length; i++) {
    if (patternName[i] == '(') {
      first = i;
      break;
    }
  }
  return patternName.slice(0, first - 1);
}

function getUserName(patternName: string): string[] {
  let first = -1;
  for (let i = 0; i < patternName.length; i++) {
    if (patternName[i] == '(') {
      first = i;
      break;
    }
  }
  return (first == -1) ? ['', ''] : patternName.slice(first + 9, patternName.length - 1).split(' ').slice(1);
}

function translateSequence(
  sequence: string,
  bases: any,
  ptoLinkages: any,
  startModification: any,
  endModification: any,
  firstPtoExist: boolean): string {
  let counter: number = -1;
  let mainSequence = sequence.replace(/[AUGC]/g, function(x: string) {
    counter++;
    const indexOfSymbol = axolabsMap['RNA']['symbols'].indexOf(x);
    const symbol = axolabsMap[bases[counter].value]['symbols'][indexOfSymbol];
    return (ptoLinkages[counter].value) ? symbol + 's' : symbol;
  });
  if (mainSequence.slice(0, 5).split('mU').length == 3)
    mainSequence = '(uu)' + mainSequence.slice(4);
  if (mainSequence.slice(mainSequence.length - 7).split('mU').length == 3)
    mainSequence = mainSequence.slice(0, mainSequence.length - 4) + '(uu)';
  return startModification.value + (firstPtoExist ? 's' : '') + mainSequence + endModification.value;
}

function addColumnWithIds(tableName: string, columnName: string, patternName: string) {
  const nameOfNewColumn = 'ID ' + patternName;
  const columns = grok.shell.table(tableName).columns;
  if (columns.contains(nameOfNewColumn))
    columns.remove(nameOfNewColumn);
  const columnWithIds = columns.byName(columnName);
  return columns.addNewString(nameOfNewColumn).init((i: number) => columnWithIds.get(i) + '_' + patternName);
}

function addColumnWithTranslatedSequences(
  tableName: string,
  columnName: string,
  bases: any,
  ptoLinkages: any,
  startModification: any,
  endModification: any,
  firstPtoExist: boolean) {
  const nameOfNewColumn = 'Axolabs ' + columnName;
  const columns = grok.shell.table(tableName).columns;
  if (columns.contains(nameOfNewColumn))
    columns.remove(nameOfNewColumn);
  const columnWithInputSequences = columns.byName(columnName);
  return columns.addNewString(nameOfNewColumn).init((i: number) => {
    return translateSequence(columnWithInputSequences.getString(i),
      bases, ptoLinkages, startModification, endModification, firstPtoExist);
  });
}

export function defineAxolabsPattern() {
  const enumerateModifications = [defaultBase];
  let maximalSsLength = defaultSequenceLength;
  let maximalAsLength = defaultSequenceLength;

  function updateAsModification() {
    asModificationItems.innerHTML = '';
    asPtoLinkages = asPtoLinkages.concat(Array(maximalAsLength - asBases.length).fill(fullyPto));
    asBases = asBases.concat(Array(maximalAsLength - asBases.length).fill(sequenceBase));
    let nucleotideCounter = 0;
    for (let i = 0; i < asLength.value; i++) {
      asPtoLinkages[i] = ui.boolInput('', asPtoLinkages[i].value, () => {
        updateSvgScheme();
        updateOutputExamples();
      });
      asBases[i] = ui.choiceInput('', asBases[i].value, baseChoices, (v: string) => {
        if (!enumerateModifications.includes(v)) {
          enumerateModifications.push(v);
          isEnumerateModificationsDiv.append(ui.boolInput(v, true, (boolV: boolean) => {
            if (boolV) {
              if (!enumerateModifications.includes(v))
                enumerateModifications.push(v);
            } else {
              const index = enumerateModifications.indexOf(v, 0);
              if (index > -1)
                enumerateModifications.splice(index, 1);
            }
            updateSvgScheme();
          }).root);
        }
        updateAsModification();
        updateSvgScheme();
        updateOutputExamples();
      });
      if (asBases[i].value.slice(-3) != '(o)')
        nucleotideCounter++;

      asModificationItems.append(
        ui.divH([
          ui.div([ui.label(asBases[i].value.slice(-3) == '(o)' ? '' : String(nucleotideCounter))],
            {style: {width: '20px'}})!,
          ui.block75([asBases[i]])!,
          ui.div([asPtoLinkages[i]])!,
        ], {style: {alignItems: 'center'}}),
      );
    }
  }

  function updateSsModification() {
    ssModificationItems.innerHTML = '';
    ssPtoLinkages = ssPtoLinkages.concat(Array(maximalSsLength - ssBases.length).fill(fullyPto));
    ssBases = ssBases.concat(Array(maximalSsLength - ssBases.length).fill(sequenceBase));
    let nucleotideCounter = 0;
    for (let i = 0; i < ssLength.value; i++) {
      ssPtoLinkages[i] = ui.boolInput('', ssPtoLinkages[i].value, () => {
        updateSvgScheme();
        updateOutputExamples();
      });
      ssBases[i] = ui.choiceInput('', ssBases[i].value, baseChoices, (v: string) => {
        if (!enumerateModifications.includes(v)) {
          enumerateModifications.push(v);
          isEnumerateModificationsDiv.append(ui.boolInput(v, true, (boolV: boolean) => {
            if (boolV) {
              if (!enumerateModifications.includes(v))
                enumerateModifications.push(v);
            } else {
              const index = enumerateModifications.indexOf(v, 0);
              if (index > -1)
                enumerateModifications.splice(index, 1);
            }
            updateSvgScheme();
          }).root);
        }
        updateSsModification();
        updateSvgScheme();
        updateOutputExamples();
      });
      if (ssBases[i].value.slice(-3) != '(o)')
        nucleotideCounter++;

      ssModificationItems.append(
        ui.divH([
          ui.div([ui.label(ssBases[i].value.slice(-3) == '(o)' ? '' : String(nucleotideCounter))],
            {style: {width: '20px'}})!,
          ui.block75([ssBases[i]])!,
          ui.div([ssPtoLinkages[i]])!,
        ], {style: {alignItems: 'center'}}),
      );
    }
  }

  function updateUiForNewSequenceLength() {
    if (ssLength.value < maximalValidSequenceLength && asLength.value < maximalValidSequenceLength) {
      if (ssLength.value > maximalSsLength)
        maximalSsLength = ssLength.value;
      if (asLength.value > maximalAsLength)
        maximalAsLength = asLength.value;
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
    for (let i = 0; i < ssPtoLinkages.length; i++)
      ssPtoLinkages[i].value = newPtoValue;

    for (let i = 0; i < asPtoLinkages.length; i++)
      asPtoLinkages[i].value = newPtoValue;

    updateSvgScheme();
  }

  function updateBases(newBasisValue: string): void {
    for (let i = 0; i < ssBases.length; i++)
      ssBases[i].value = newBasisValue;

    for (let i = 0; i < asBases.length; i++)
      asBases[i].value = newBasisValue;

    updateSvgScheme();
  }

  function updateInputExamples() {
    ssInputExample.value = generateExample(ssLength.value, sequenceBase.value);
    if (createAsStrand.value)
      asInputExample.value = generateExample(asLength.value, sequenceBase.value);
  }

  function updateOutputExamples() {
    ssOutputExample.value = translateSequence(
      ssInputExample.value, ssBases, ssPtoLinkages, ssFiveModification, ssThreeModification, firstSsPto.value);
    if (createAsStrand.value) {
      asOutputExample.value = translateSequence(
        asInputExample.value, asBases, asPtoLinkages, asFiveModification, asThreeModification, firstAsPto.value);
    }
  }

  function updateSvgScheme() {
    svgDiv.innerHTML = '';
    svgDiv.append(
      ui.span([
        drawAxolabsPattern(
          getShortName(saveAs.value),
          createAsStrand.value,
          ssBases.slice(0, ssLength.value).map((e) => e.value),
          asBases.slice(0, asLength.value).map((e) => e.value),
          [firstSsPto.value].concat(ssPtoLinkages.slice(0, ssLength.value).map((e) => e.value)),
          [firstAsPto.value].concat(asPtoLinkages.slice(0, asLength.value).map((e) => e.value)),
          ssThreeModification.value,
          ssFiveModification.value,
          asThreeModification.value,
          asFiveModification.value,
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
      sequenceBase.value = detectDefaultBasis(obj['asBases'].concat(obj['ssBases']));
      createAsStrand.value = (obj['asBases'].length > 0);
      saveAs.value = newName;

      ssBases = [];
      for (let i = 0; i < obj['ssBases'].length; i++)
        ssBases.push(ui.choiceInput('', obj['ssBases'][i], baseChoices));

      asBases = [];
      for (let i = 0; i < obj['asBases'].length; i++)
        asBases.push(ui.choiceInput('', obj['asBases'][i], baseChoices));

      firstSsPto.value = obj['ssPtoLinkages'][0];
      ssPtoLinkages = [];
      for (let i = 1; i < obj['ssPtoLinkages'].length; i++)
        ssPtoLinkages.push(ui.boolInput('', obj['ssPtoLinkages'][i]));

      firstAsPto.value = obj['asPtoLinkages'][0];
      asPtoLinkages = [];
      for (let i = 1; i < obj['asPtoLinkages'].length; i++)
        asPtoLinkages.push(ui.boolInput('', obj['asPtoLinkages'][i]));

      ssLength.value = obj['ssBases'].length;
      asLength.value = obj['asBases'].length;

      ssThreeModification.value = obj['ssThreeModification'];
      ssFiveModification.value = obj['ssFiveModification'];
      asThreeModification.value = obj['asThreeModification'];
      asFiveModification.value = obj['asFiveModification'];
      comment.value = obj['comment'];
    });
    pi.close();
  }

  function checkWhetherAllValuesInColumnHaveTheSameLength(colName: string): boolean {
    const col = tables.value.columns.byName(colName);
    let allLengthsAreTheSame = true;
    for (let i = 1; i < col.length; i++) {
      if (col.get(i - 1).length != col.get(i).length) {
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
          tables.value.columns.addNewInt('Sequences lengths in ' + colName).init((j: number) => col.get(j).length);
          grok.shell.info('Column with lengths added to \'' + tables.value.name + '\'');
          dialog.close();
          grok.shell.v = grok.shell.getTableView(tables.value.name);
        })
        .show();
    }
    return allLengthsAreTheSame;
  }

  async function getCurrentUserName(): Promise<string> {
    return await grok.dapi.users.current().then((user) => {
      return ' (created by ' + user.firstName + ' ' + user.lastName + ')';
    });
  }

  async function postPatternToUserStorage() {
    const author = await getCurrentUserName();
    if (!saveAs.stringValue.includes('(created by '))
      saveAs.value = saveAs.stringValue + author;
    return grok.dapi.userDataStorage.postValue(
      userStorageKey,
      saveAs.stringValue,
      JSON.stringify({
        'ssBases': ssBases.slice(0, ssLength.value).map((e) => e.value),
        'asBases': asBases.slice(0, asLength.value).map((e) => e.value),
        'ssPtoLinkages': [firstSsPto.value].concat(ssPtoLinkages.slice(0, ssLength.value).map((e) => e.value)),
        'asPtoLinkages': [firstAsPto.value].concat(asPtoLinkages.slice(0, asLength.value).map((e) => e.value)),
        'ssThreeModification': ssThreeModification.value,
        'ssFiveModification': ssFiveModification.value,
        'asThreeModification': asThreeModification.value,
        'asFiveModification': asFiveModification.value,
        'comment': comment.value,
      }),
      false,
    ).then(() => grok.shell.info('Pattern \'' + saveAs.value + '\' was successfully uploaded!'));
  }

  async function updatePatternsList() {
    grok.dapi.userDataStorage.get(userStorageKey, false).then(async (entities) => {
      const lstMy: string[] = [];
      const lstOthers: string[] = [];

      for (const ent of Object.keys(entities)) {
        if (await isCurrentUserCreatedThisPattern(ent))
          lstOthers.push(ent);
        else
          lstMy.push(getShortName(ent));
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

  const inputSsColumnDiv = ui.div([]);
  const inputAsColumnDiv = ui.div([]);
  const inputIdColumnDiv = ui.div([]);
  const ssModificationItems = ui.div([]);
  const asModificationItems = ui.div([]);
  const svgDiv = ui.div([]);
  const asExampleDiv = ui.div([]);
  const appAxolabsDescription = ui.div([]);
  const loadPatternDiv = ui.div([]);
  const asModificationDiv = ui.div([]);
  const firstAsPtoDiv = ui.div([]);
  const isEnumerateModificationsDiv = ui.divH([ui.boolInput(defaultBase, true, (v: boolean) => {
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
  }).root]);

  let ssBases = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices));
  let asBases = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices));
  let ssPtoLinkages = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto));
  let asPtoLinkages = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto));

  const ssLength = ui.intInput('SS Length', defaultSequenceLength, () => updateUiForNewSequenceLength());
  const asLength = ui.intInput('AS Length', defaultSequenceLength, () => updateUiForNewSequenceLength());
  const asLengthDiv = ui.div([asLength.root]);

  function validateSsColumn(colName: string) {
    const allLengthsAreTheSame: boolean = checkWhetherAllValuesInColumnHaveTheSameLength(colName);
    const firstSequence = tables.value.columns.byName(colName).get(0);
    if (allLengthsAreTheSame && firstSequence.length != ssLength.value)
      ssLength.value = tables.value.columns.byName(colName).get(0).length;
    ssInputExample.value = firstSequence;
  }

  function validateAsColumn(colName: string) {
    const allLengthsAreTheSame: boolean = checkWhetherAllValuesInColumnHaveTheSameLength(colName);
    const firstSequence = tables.value.columns.byName(colName).get(0);
    if (allLengthsAreTheSame && firstSequence.length != asLength.value)
      asLength.value = tables.value.columns.byName(colName).get(0).length;
    asLengthDiv.innerHTML = '';
    asLengthDiv.append(asLength.root);
    asInputExample.value = firstSequence;
  }

  function validateIdsColumn(colName: string) {
    const col = tables.value.columns.byName(colName);
    if (col.type != DG.TYPE.INT)
      grok.shell.error('Column should contain integers only');
    else if (col.categories.length < col.length) {
      const duplicates = findDuplicates(col.getRawData());
      ui.dialog('Non-unique IDs')
        .add(ui.divText('Press \'OK\' to select rows with non-unique values'))
        .onOK(() => {
          const selection = tables.value.selection;
          selection.init((i: number) => duplicates.indexOf(col.get(i)) > -1);
          grok.shell.v = grok.shell.getTableView(tables.value.name);
          grok.shell.info('Rows are selected in table \'' + tables.value.name + '\'');
        })
        .show();
    }
  }

  const tables = ui.tableInput('Tables', grok.shell.tables[0], grok.shell.tables, (t: DG.DataFrame) => {
    inputSsColumn = ui.choiceInput('SS Column', '', t.columns.names(), (colName: string) => validateSsColumn(colName));
    inputSsColumnDiv.innerHTML = '';
    inputSsColumnDiv.append(inputSsColumn.root);
    inputAsColumn = ui.choiceInput('AS Column', '', t.columns.names(), (colName: string) => validateAsColumn(colName));
    inputAsColumnDiv.innerHTML = '';
    inputAsColumnDiv.append(inputAsColumn.root);
    inputIdColumn = ui.choiceInput('ID Column', '', t.columns.names(), (colName: string) => validateIdsColumn(colName));
    inputIdColumnDiv.innerHTML = '';
    inputIdColumnDiv.append(inputIdColumn.root);
  });

  let inputSsColumn = ui.choiceInput('SS Column', '', []);
  inputSsColumnDiv.append(inputSsColumn.root);
  let inputAsColumn = ui.choiceInput('AS Column', '', []);
  inputAsColumnDiv.append(inputAsColumn.root);
  let inputIdColumn = ui.choiceInput('ID Column', '', []);
  inputIdColumnDiv.append(inputIdColumn.root);

  updatePatternsList();

  const sequenceBase = ui.choiceInput('Sequence Basis', defaultBase, baseChoices, (v: string) => {
    updateBases(v);
    updateOutputExamples();
  });

  const fullyPto = ui.boolInput('Fully PTO', defaultPto, (v: boolean) => {
    firstSsPto.value = v;
    firstAsPto.value = v;
    updatePto(v);
    updateOutputExamples();
  });

  const firstSsPto = ui.boolInput('First SS PTO', fullyPto.value, () => updateSvgScheme());
  const firstAsPto = ui.boolInput('First AS PTO', fullyPto.value, () => updateSvgScheme());
  firstAsPtoDiv.append(firstAsPto.root);

  const createAsStrand = ui.boolInput('Create AS Strand', true, (v: boolean) => {
    asModificationSection.hidden = (!v);
    inputAsColumnDiv.hidden = (!v);
    asLengthDiv.hidden = (!v);
    asModificationDiv.hidden = (!v);
    asExampleDiv.hidden = (!v);
    firstAsPtoDiv.hidden = (!v);
    updateSvgScheme();
  });

  const saveAs = ui.textInput('Save As', 'Pattern Name', () => updateSvgScheme());
  saveAs.setTooltip('Name Of New Pattern');

  const ssThreeModification = ui.stringInput('SS 3\' Modification', '', () => {
    updateSvgScheme();
    updateOutputExamples();
  });
  ssThreeModification.setTooltip('Additional SS 3\' Modification');

  const ssFiveModification = ui.stringInput('SS 5\' Modification', '', () => {
    updateSvgScheme();
    updateOutputExamples();
  });
  ssFiveModification.setTooltip('Additional SS 5\' Modification');

  const asThreeModification = ui.stringInput('AS 3\' Modification', '', () => {
    updateSvgScheme();
    updateOutputExamples();
  });
  asThreeModification.setTooltip('Additional AS 3\' Modification');

  const asFiveModification = ui.stringInput('AS 5\' Modification', '', () => {
    updateSvgScheme();
    updateOutputExamples();
  });
  asFiveModification.setTooltip('Additional AS 5\' Modification');

  asModificationDiv.append(asThreeModification.root);
  asModificationDiv.append(asFiveModification.root);

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
    if (inputSsColumn.value == null || (createAsStrand.value && inputAsColumn.value == null))
      grok.shell.info('Please select table and columns on which to apply pattern');
    else if (ssLength.value != ssInputExample.value.length || asLength.value != asInputExample.value.length) {
      const dialog = ui.dialog('Length Mismatch');
      $(dialog.getButton('OK')).hide();
      dialog
        .add(ui.divText('Length of sequences in columns doesn\'t match entered length. Update length value?'))
        .addButton('YES', () => {
          ssLength.value = tables.value.columns.byName(inputSsColumn.value).getString(0).length;
          asLength.value = tables.value.columns.byName(inputAsColumn.value).getString(0).length;
          dialog.close();
        })
        .show();
    } else {
      if (inputIdColumn.value != null)
        addColumnWithIds(tables.value.name, inputIdColumn.value, getShortName(saveAs.value));
      addColumnWithTranslatedSequences(
        tables.value.name, inputSsColumn.value, ssBases, ssPtoLinkages,
        ssFiveModification, ssThreeModification, firstSsPto.value);
      if (createAsStrand.value) {
        addColumnWithTranslatedSequences(
          tables.value.name, inputAsColumn.value, asBases, asPtoLinkages,
          asFiveModification, asThreeModification, firstAsPto.value);
      }
      grok.shell.v = grok.shell.getTableView(tables.value.name);
      grok.shell.info(((createAsStrand.value) ? 'Columns were' : 'Column was') +
      ' added to table \'' + tables.value.name + '\'');
    }
  });

  const ssInputExample = ui.textInput('Sense Strand', generateExample(ssLength.value, sequenceBase.value), () => {
    ssOutputExample.value = translateSequence(ssInputExample.value, ssBases, ssPtoLinkages,
      ssFiveModification, ssThreeModification, firstSsPto.value);
  });
  const ssOutputExample = ui.textInput(' ', translateSequence(
    ssInputExample.value, ssBases, ssPtoLinkages, ssThreeModification, ssFiveModification, firstSsPto.value));
  (ssInputExample.input as HTMLElement).style.resize = 'none';
  (ssInputExample.input as HTMLElement).style.minWidth = exampleMinWidth;
  (ssOutputExample.input as HTMLElement).style.resize = 'none';
  (ssOutputExample.input as HTMLElement).style.minWidth = exampleMinWidth;
  // @ts-ignore
  ssOutputExample.input.disabled = 'true';
  ssOutputExample.root.append(
    ui.div([
      ui.button(ui.iconFA('copy', () => {}), () => {
        navigator.clipboard.writeText(ssOutputExample.value).then(() =>
          grok.shell.info('Sequence was copied to clipboard'));
      }),
    ], 'ui-input-options'),
  );

  const asInputExample = ui.textInput('Antisense Strand', generateExample(asLength.value, sequenceBase.value), () => {
    asOutputExample.value = translateSequence(
      asInputExample.value, asBases, asPtoLinkages, asFiveModification, asThreeModification, firstSsPto.value);
  });
  const asOutputExample = ui.textInput(' ', translateSequence(
    asInputExample.value, asBases, asPtoLinkages, asFiveModification, asThreeModification, firstSsPto.value));
  (asInputExample.input as HTMLElement).style.resize = 'none';
  (asInputExample.input as HTMLElement).style.minWidth = exampleMinWidth;
  (asOutputExample.input as HTMLElement).style.resize = 'none';
  (asOutputExample.input as HTMLElement).style.minWidth = exampleMinWidth;
  // @ts-ignore
  asOutputExample.input.disabled = 'true';
  asOutputExample.root.append(
    ui.div([
      ui.button(ui.iconFA('copy', () => {}), () => {
        navigator.clipboard.writeText(asOutputExample.value).then(() =>
          grok.shell.info('Sequence was copied to clipboard'));
      }),
    ], 'ui-input-options'),
  );
  asExampleDiv.append(asInputExample.root);
  asExampleDiv.append(asOutputExample.root);

  updateUiForNewSequenceLength();

  const exampleSection = ui.div([
    ui.h1('Example'),
    ssInputExample.root,
    ssOutputExample.root,
    asExampleDiv,
  ], 'ui-form');

  const inputsSection = ui.div([
    ui.h1('Inputs'),
    ui.divH([
      tables.root,
      inputSsColumnDiv,
    ]),
    ui.divH([
      inputAsColumnDiv,
      inputIdColumnDiv,
    ]),
    ui.buttonsInput([
      convertSequenceButton,
    ]),
  ], 'ui-form');

  const mainSection = ui.panel([
    ui.block([
      svgDiv,
    ], {style: {overflowX: 'scroll'}}),
    ui.button('Download', () => svg.saveSvgAsPng(document.getElementById('mySvg'), saveAs.value)),
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
            ssLength.root,
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
            firstSsPto.root,
            firstAsPtoDiv,
            ssFiveModification.root,
            ssThreeModification.root,
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
    ssModificationItems,
  ])!;

  const asModificationSection = ui.panel([
    ui.h1('Antisense Strand'),
    ui.divH([
      ui.div([ui.divText('#')], {style: {width: '20px'}})!,
      ui.block75([ui.divText('Modification')])!,
      ui.div([ui.divText('PTO')], {style: {paddingRight: '8px'}})!,
    ]),
    asModificationItems,
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
