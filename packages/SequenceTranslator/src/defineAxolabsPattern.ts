/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// @ts-ignore
import * as svg from 'save-svg-as-png';

import {drawAxolabsPattern} from "./drawAxolabsPattern";
import {axolabsMap} from "./axolabsMap";
/*
SS - sense strand of DNA;
AS - antisense strand of DNA;
base - indicates how to translate input nucleotide;
PTO - indicates whether oligonucleotide is phosphorothioated (ps linkage);
SS/AS Modification - sections of dialog for changing base and PTO statuses of one nucleotide at once: for SS or for AS;
*/

const baseChoices: string[] = Object.keys(axolabsMap);
const defaultBase: string = baseChoices[0];
const defaultPto: boolean = false;
const defaultSequenceLength: number = 23;
const maximalValidSequenceLength: number = 35;
const userStorageKey: string = 'SequenceTranslator';
const exampleMinWidth: string = '400px';
const modificationSectionMaxWidth = '220px';

function translateSequence(sequence: string, bases: any, ptoLinkages: any, modification: any) {
  let counter: number = -1;
  return sequence.replace(/[AUGC]/g, function (x: string) {
    counter++;
    let indexOfSymbol = axolabsMap["RNA"]["symbols"].indexOf(x);
    let symbol = axolabsMap[bases[counter].value]["symbols"][indexOfSymbol];
    return (ptoLinkages[counter].value) ? symbol + 's' : symbol;
  }) + modification.value;
}

function addColumnWithTranslatedSequences(tableName: string, columnName: string, bases: any, ptoLinkages: any, modification: any) {
  const nameOfNewColumn = 'Axolabs ' + columnName;
  let columns = grok.shell.table(tableName).columns;
  if (columns.contains(nameOfNewColumn))
    columns.remove(nameOfNewColumn);
  const columnWithInputSequences = columns.byName(columnName);
  return columns.addNewString(nameOfNewColumn).init((i: number) => {
    return translateSequence(columnWithInputSequences.getString(i), bases, ptoLinkages, modification);
  });
}

export function defineAxolabsPattern() {

  function updateAsModification() {
    asModificationItems.innerHTML = '';
    asPtoLinkages = (asLength.value > asPtoLinkages.length) ?
      asPtoLinkages.concat(Array(asLength.value - asPtoLinkages.length).fill(fullyPto)) :
      asPtoLinkages.slice(asPtoLinkages.length - asLength.value);
    asBases = (asLength.value > asBases.length) ?
      asBases.concat(Array(asLength.value - asBases.length).fill(sequenceBase)) :
      asBases.slice(asBases.length - asLength.value);
    for (let i = 0; i < asLength.value; i++) {
      asPtoLinkages[i] = ui.boolInput('', asPtoLinkages[i].value, () => {
        updateSvgScheme();
        updateExamples();
      });
      asBases[i] = ui.choiceInput('', asBases[i].value, baseChoices, () => {
        updateSvgScheme();
        updateExamples();
      });
      asModificationItems.append(
        ui.divH([
          ui.block25([ui.label((i + 1).toString())])!,
          ui.block50([asBases[i]])!,
          ui.block25([asPtoLinkages[i]])!
        ], {style: {alignItems: "center"}})
      );
    }
  }

  function updateSsModification() {
    ssModificationItems.innerHTML = '';
    ssPtoLinkages = (ssLength.value > ssPtoLinkages.length) ?
      ssPtoLinkages.concat(Array(ssLength.value - ssPtoLinkages.length).fill(fullyPto)) :
      ssPtoLinkages.slice(ssPtoLinkages.length - ssLength.value);
    ssBases = (ssLength.value > ssBases.length) ?
      ssBases.concat(Array(ssLength.value - ssBases.length).fill(sequenceBase)) :
      ssBases.slice(ssBases.length - ssLength.value);
    for (let i = 0; i < ssLength.value; i++) {
      ssPtoLinkages[i] = ui.boolInput('', ssPtoLinkages[i].value, () => {
        updateSvgScheme();
        updateExamples();
      });
      ssBases[i] = ui.choiceInput('', ssBases[i].value, baseChoices, () => {
        updateSvgScheme();
        updateExamples();
      });
      ssModificationItems.append(
        ui.divH([
          ui.block25([ui.label((i + 1).toString())])!,
          ui.block50([ssBases[i]])!,
          ui.block25([ssPtoLinkages[i]])!
        ], {style: {alignItems: "center"}})
      );
    }
  }

  function updateUiForNewSequenceLength() {
    if (ssLength.value < maximalValidSequenceLength && asLength.value < maximalValidSequenceLength) {
      updateSsModification();
      updateAsModification();
      updateSvgScheme();
    } else {
      ui.dialog('Sequence length is out of range')
        .add(ui.divText('Sequence length should be less than ' + maximalValidSequenceLength.toString() + ' due to UI constrains.'))
        .add(ui.divText('Please change sequence length in order to define new pattern.'))
        .show();
    }
  }

  function updatePto(newPtoValue: boolean) {
    for (let i = 0; i < ssPtoLinkages.length; i++) {
      ssPtoLinkages[i].value = newPtoValue;
    }
    for (let i = 0; i < asPtoLinkages.length; i++) {
      asPtoLinkages[i].value = newPtoValue;
    }
    updateSvgScheme();
  }

  function updateBases(newBasisValue: string) {
    for (let i = 0; i < ssBases.length; i++) {
      ssBases[i].value = newBasisValue;
    }
    for (let i = 0; i < asBases.length; i++) {
      asBases[i].value = newBasisValue;
    }
    updateSvgScheme();
  }

  function updateExamples() {
    ssResult.value = translateSequence(ssInput.value, ssBases, ssPtoLinkages, threeModification);
    if (createAsStrand.value)
      asResult.value = translateSequence(asInput.value, asBases, asPtoLinkages, fiveModification);
  }

  function updateSvgScheme() {
    svgDiv.innerHTML = '';
    svgDiv.append(
      ui.span([
        drawAxolabsPattern(
          saveAs.value,
          createAsStrand.value,
          ssBases.slice(0, ssLength.value).map((e) => e.value),
          asBases.slice(0, asLength.value).map((e) => e.value),
          ssPtoLinkages.slice(0, ssLength.value).map((e) => e.value),
          asPtoLinkages.slice(0, asLength.value).map((e) => e.value),
          threeModification.value,
          fiveModification.value
        )
      ])
    );
  }

  function parsePatternAndUpdateUi(newName: string) {
    grok.dapi.userDataStorage.get(userStorageKey, false).then((entities) => {
      // @ts-ignore
      let obj = JSON.parse(entities[newName]);
      ssLength.value = obj['ssBases'].length;
      asLength.value = obj['asBases'].length;
      createAsStrand.value = (asLength.value > 0);
      saveAs.value = newName;

      ssBases = [];
      for (let i = 0; i < obj['ssBases'].length; i++)
        ssBases.push(ui.choiceInput('', obj['ssBases'][i], baseChoices));

      asBases = [];
      for (let i = 0; i < obj['asBases'].length; i++)
        asBases.push(ui.choiceInput('', obj['asBases'][i], baseChoices));

      ssPtoLinkages = [];
      for (let i = 0; i < obj['ssPtoLinkages'].length; i++)
        ssPtoLinkages.push(ui.boolInput('', obj['ssPtoLinkages'][i]));

      asPtoLinkages = [];
      for (let i = 0; i < obj['asPtoLinkages'].length; i++)
        asPtoLinkages.push(ui.boolInput('', obj['asPtoLinkages'][i]));

      threeModification.value = obj['threeModification'];
      fiveModification.value = obj['fiveModification'];

      updateSvgScheme();
      updateAsModification();
      updateSsModification();
    });
  }

  function checkWhetherAllValuesInColumnHaveTheSameLength(colName: string): boolean {
    let col = grok.shell.table(tables.value).columns.byName(colName);
    let allLengthsAreTheSame = true;
    for (let i = 1; i < col.length; i++) {
      if (col.get(i - 1).length != col.get(i).length) {
        allLengthsAreTheSame = false;
        break;
      }
    }
    if (!allLengthsAreTheSame) {
      let dialog = ui.dialog('Sequences lengths mismatch');
      $(dialog.getButton('OK')).hide();
      dialog
        .add(ui.divText('The sequence length should match the number of Raw sequences in the input file'))
        .add(ui.divText("'ADD COLUMN' to see sequences lengths"))
        .addButton('ADD COLUMN', () => {
          grok.shell.table(tables.value).columns.addNewInt('Sequences lengths in ' + colName).init((j) => col.get(j).length);
          grok.shell.info("Column with lengths added to '" + tables.value + "'");
          dialog.close();
        })
        .show();
    }
    return allLengthsAreTheSame;
  }

  function postPatternToUserStorage() {
    grok.dapi.userDataStorage.postValue(
      userStorageKey,
      saveAs.value,
      JSON.stringify({
        "ssBases": ssBases.slice(0, ssLength.value).map((e) => e.value),
        "asBases": asBases.slice(0, asLength.value).map((e) => e.value),
        "ssPtoLinkages": ssPtoLinkages.slice(0, ssLength.value).map((e) => e.value),
        "asPtoLinkages": asPtoLinkages.slice(0, asLength.value).map((e) => e.value),
        "threeModification": threeModification.value,
        "fiveModification": fiveModification.value
      }),
      false
    ).then(() => grok.shell.info('Pattern ' + saveAs.value + ' was successfully uploaded!'));
  }

  function updatePatternsList() {
    grok.dapi.userDataStorage.get(userStorageKey, false).then((entities) => {
      let loadPattern = ui.choiceInput('Load Pattern', '', Object.keys(entities), (v: string) => {
        parsePatternAndUpdateUi(v)
      });
      loadPattern.setTooltip('Apply Existing Pattern');
      loadPatternDiv.innerHTML = '';
      loadPatternDiv.append(
        loadPattern.root
      );
      loadPattern.root.append(
        ui.button(ui.iconFA('trash-alt', () => {
        }), () => {
          grok.dapi.userDataStorage.remove(userStorageKey, loadPattern.value, false)
            .then(() => grok.shell.info("Pattern '" + loadPattern.value + "' deleted"));
          updatePatternsList();
        })
      );
    });
  }

  function savePattern() {
    grok.dapi.userDataStorage.get(userStorageKey, false).then((entities) => {
      if (Object.keys(entities).includes(saveAs.value)) {
        let dialog = ui.dialog('Pattern already exists');
        $(dialog.getButton('OK')).hide();
        dialog
          .add(ui.divText("Pattern name '" + saveAs.value + "' already exists."))
          .add(ui.divText('Replace pattern?'))
          .addButton('YES', () => {
            grok.dapi.userDataStorage.remove(userStorageKey, saveAs.value, false)
              .then(() => postPatternToUserStorage());
            dialog.close();
          })
          .show();
      } else {
        postPatternToUserStorage();
      }
    });
    updatePatternsList();
  }

  let inputSsColumnDiv = ui.div([]),
    inputAsColumnDiv = ui.div([]),
    ssModificationItems = ui.div([]),
    asModificationItems = ui.div([]),
    svgDiv = ui.div([]),
    asExampleDiv = ui.div([]),
    appAxolabsDescription = ui.div([]),
    loadPatternDiv = ui.div([]);

  let ssBases = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices)),
    asBases = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices)),
    ssPtoLinkages = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto)),
    asPtoLinkages = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto));

  let ssLength = ui.intInput('SS Length', defaultSequenceLength, () => updateUiForNewSequenceLength());
  let asLength = ui.intInput('AS Length', defaultSequenceLength, () => updateUiForNewSequenceLength());
  let asLengthDiv = ui.div([asLength.root]);

  function f1(n: string) {
    let allLengthsAreTheSame: boolean = checkWhetherAllValuesInColumnHaveTheSameLength(n);
    const firstSequence = grok.shell.table(tables.value).columns.byName(n).get(0);
    if (allLengthsAreTheSame && firstSequence.length != ssLength.value)
      ssLength.value = grok.shell.table(tables.value).columns.byName(n).get(0).length;
      ssInput.value = firstSequence;
  }

  function f2(n: string) {
    let allLengthsAreTheSame: boolean = checkWhetherAllValuesInColumnHaveTheSameLength(n);
    const firstSequence = grok.shell.table(tables.value).columns.byName(n).get(0);
    if (allLengthsAreTheSame && firstSequence.length != asLength.value)
      asLength.value = grok.shell.table(tables.value).columns.byName(n).get(0).length;
      asLengthDiv.innerHTML = '';
      asLengthDiv.append(asLength.root);
      asInput.value = firstSequence;
  }

  let tables = ui.choiceInput('Tables', '', grok.shell.tableNames, () => {
    inputSsColumn = ui.choiceInput('SS Column', '', grok.shell.table(tables.value).columns.names(), (n: string) => f1(n));
    inputSsColumnDiv.innerHTML = '';
    inputSsColumnDiv.append(inputSsColumn.root);
    inputAsColumn = ui.choiceInput('AS Column', '', grok.shell.table(tables.value).columns.names(), (n: string) => f2(n));
    inputAsColumnDiv.innerHTML = '';
    inputAsColumnDiv.append(inputAsColumn.root);
  });

  let inputSsColumn = ui.choiceInput('SS Column', '', []);
  inputSsColumnDiv.append(inputSsColumn.root);
  let inputAsColumn = ui.choiceInput('AS Column', '', []);
  inputAsColumnDiv.append(inputAsColumn.root);

  updatePatternsList();

  let sequenceBase = ui.choiceInput('Sequence Basis', defaultBase, baseChoices, (v: string) => {
    updateBases(v);
    updateExamples();
  });

  let fullyPto = ui.boolInput('Fully PTO', defaultPto, (v: boolean) => {
    updatePto(v);
    updateExamples();
  });

  let createAsStrand = ui.boolInput('Create AS Strand', true, (v: boolean) => {
    asModificationSection.hidden = (!v);
    inputAsColumnDiv.hidden = (!v);
    asLengthDiv.hidden = (!v);
    asExampleDiv.hidden = (!v);
    updateSvgScheme();
  });

  let saveAs = ui.stringInput('Save As', '', () => updateSvgScheme());
  saveAs.setTooltip('Name Of New Pattern');

  let threeModification = ui.stringInput("Additional 3' Modification", "", () => {
    updateSvgScheme();
    updateExamples();
  });
  threeModification.setTooltip("Additional 3' Modification");
  let fiveModification = ui.stringInput("Additional 5' Modification", "", () => {
    updateSvgScheme();
    updateExamples();
  });
  fiveModification.setTooltip("Additional 5' Modification");

  updateUiForNewSequenceLength();

  let savePatternButton = ui.button('Save', () => {
    if (saveAs.value != '') {
      savePattern();
    } else {
      let name = ui.stringInput('Enter name', '');
      ui.dialog('Pattern name')
        .add(name.root)
        .onOK(() => {
          saveAs.value = name.value;
          savePattern();
        })
        .show();
    }
  });

  let convertSequenceButton = ui.button('Convert Sequences', () => {
    if (inputSsColumn.value == null || (createAsStrand.value && inputAsColumn.value == null))
      grok.shell.info("Please select table and columns on which to apply pattern");
    else {
      addColumnWithTranslatedSequences(tables.value, inputSsColumn.value, ssBases, ssPtoLinkages, threeModification);
      if (createAsStrand.value)
        addColumnWithTranslatedSequences(tables.value, inputAsColumn.value, asBases, asPtoLinkages, fiveModification);
      grok.shell.info(((createAsStrand.value) ? "Columns were" : "Column was") + " added to table '" + tables.value + "'");
    }
  });

  let ssInput = ui.stringInput('SS', '',() => ssResult.value = translateSequence(ssInput.value, ssBases, ssPtoLinkages, threeModification));
  let ssResult = ui.stringInput(' ', '');
  // @ts-ignore
  ssInput.input.style.resize = 'none';
  // @ts-ignore
  ssInput.input.style.minWidth = exampleMinWidth;
  // @ts-ignore
  ssResult.input.style.resize = 'none';
  // @ts-ignore
  ssResult.input.style.minWidth = exampleMinWidth;
  // @ts-ignore
  ssResult.input.disabled = 'true';
  ssResult.root.append(
    ui.button(ui.iconFA('copy', () => {}), () => {
      navigator.clipboard.writeText(ssResult.value).then(() => grok.shell.info('Sequence was copied to clipboard'));
    })
  );

  let asInput = ui.stringInput('AS','',() => asResult.value = translateSequence(asInput.value, asBases, asPtoLinkages, fiveModification));
  let asResult = ui.stringInput(' ','');
  // @ts-ignore
  asInput.input.style.resize = 'none';
  // @ts-ignore
  asInput.input.style.minWidth = exampleMinWidth;
  // @ts-ignore
  asResult.input.style.resize = 'none';
  // @ts-ignore
  asResult.input.style.minWidth = exampleMinWidth;
  // @ts-ignore
  asResult.input.disabled = 'true';
  asResult.root.append(
    ui.button(ui.iconFA('copy', () => {}), () => {
      navigator.clipboard.writeText(asResult.value).then(() => grok.shell.info('Sequence was copied to clipboard'));
    })
  );
  asExampleDiv.append(asInput.root);
  asExampleDiv.append(asResult.root);

  let patternDesignSection = ui.panel([
    ui.divH([
      ui.div([
        ui.divH([
          ui.h1('Pattern'),
          ui.iconFA('question-circle',() => {
            appAxolabsDescription.innerHTML = '';
            appAxolabsDescription.append(info);
          })
        ]),
        ssLength.root,
        asLengthDiv,
        sequenceBase.root,
        fullyPto.root,
        createAsStrand.root,
        threeModification.root,
        fiveModification.root,
        loadPatternDiv,
        saveAs.root,
        ui.buttonsInput([
          savePatternButton
        ])
      ], 'ui-form'),
      ui.div([
        ui.div([
          ui.h1('Inputs'),
          tables.root,
          inputSsColumnDiv,
          inputAsColumnDiv,
          ui.buttonsInput([
            convertSequenceButton
          ])
        ], 'ui-form'),
        ui.div([
          ui.h1('Example'),
          ssInput.root,
          ssResult.root,
          asExampleDiv
        ], 'ui-form')
      ])
    ], {style: {flexWrap: 'wrap'}}),
    ui.block([
      svgDiv
    ], {style: {overflowX: 'scroll'}}),
    ui.button('Download', () => svg.saveSvgAsPng(document.getElementById('mySvg'), saveAs.value))
  ]);

  let ssModificationSection = ui.box(
    ui.panel([
      ui.h1('Sense Strand'),
      ui.divH([
        ui.block25([ui.divText('#')])!,
        ui.block50([ui.divText('Modification')])!,
        ui.block25([ui.divText('PTO')])!
      ]),
      ssModificationItems
    ])!, {style: {maxWidth: modificationSectionMaxWidth}});

  let asModificationSection = ui.box(
    ui.panel([
      ui.h1('Antisense Strand'),
      ui.divH([
        ui.block25([ui.divText('#')])!,
        ui.block50([ui.divText('Modification')])!,
        ui.block25([ui.divText('PTO')])!
      ]),
      asModificationItems
    ])!, {style: {maxWidth: modificationSectionMaxWidth}});

  let info = ui.info(
    [
      ui.divText("\n How to define new pattern:",{style: {'font-weight': 'bolder'}}),
      ui.divText("1. Choose table and columns with sense and antisense strands"),
      ui.divText("2. Choose lengths of both strands by editing checkboxes below"),
      ui.divText("3. Choose basis and PTO status for each nucleotide"),
      ui.divText("4. Set additional modifications for sequence edges"),
      ui.divText("5. Press 'Convert Sequences' button"),
      ui.divText("This will add the result column(s) to the right of the table"),
    ], 'Create and apply Axolabs translation patterns.'
  );

  return ui.splitH([
    ui.div([
      appAxolabsDescription,
      patternDesignSection!
    ])!,
    ssModificationSection,
    asModificationSection
  ]);
}