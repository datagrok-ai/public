/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {drawAxolabsPattern} from "./drawAxolabsPattern";
import {axolabsMap} from "./axolabsMap";
/*
SS - sense strand of DNA;
AS - antisense strand of DNA;
base - indicates how to translate input nucleotide;
PTO - indicates whether oligonucleotide is phosphorothioated (ps linkage);
Pattern design - section of dialog for changing length of AS and SS by editing corresponding checkboxes;
SS/AS Modification - sections of dialog for changing base and PTO statuses of one nucleotide at once: for SS or for AS;
*/

const baseChoices: string[] = Object.keys(axolabsMap);
const defaultBase: string = baseChoices[0];
const defaultPto: boolean = false;
const defaultAvailability: boolean = true;
const defaultSequenceLength: number = 23;
const maximalValidSequenceLength: number = 35;


export function defineAxolabsPattern() {

  function updateSsPattern() {
    ssPattern.innerHTML = '';
    ssAvailabilityStatuses = (sequenceLength.value > ssAvailabilityStatuses.length) ?
      ssAvailabilityStatuses.concat(Array(sequenceLength.value - ssAvailabilityStatuses.length).fill(ui.boolInput('', true))) :
      ssAvailabilityStatuses.slice(ssAvailabilityStatuses.length - sequenceLength.value);
    const ss5 = ui.divText("SS-5'");
    ss5.style.marginTop = '26px';
    ssPattern.append(ss5);
    for (let i = 0; i < ssAvailabilityStatuses.length; i++) {
      ssAvailabilityStatuses[i] = ui.boolInput('', ssAvailabilityStatuses[i].value, (v: boolean) => updateSsAvailability(i, v));
      let index = ui.divText((i + 1).toString());
      index.style.marginLeft = '7px';
      ssPattern.append(
        ui.divV([
          index,
          ssAvailabilityStatuses[i]
        ])
      );
    }
    const three = ui.divText("3'");
    three.style.marginTop = '26px';
    ssPattern.append(three);
  }

  function updateAsPattern() {
    asPattern.innerHTML = '';
    asAvailabilityStatuses = (sequenceLength.value > asAvailabilityStatuses.length) ?
      asAvailabilityStatuses.concat(Array(sequenceLength.value - asAvailabilityStatuses.length).fill(ui.boolInput('', true))) :
      asAvailabilityStatuses.slice(asAvailabilityStatuses.length - sequenceLength.value);
    const as3 = ui.divText("AS-3'");
    as3.style.marginTop = '10px';
    asPattern.append(as3);
    for (let i = 0; i < asAvailabilityStatuses.length; i++) {
      asAvailabilityStatuses[i] = ui.boolInput('', asAvailabilityStatuses[i].value, (v: boolean) => updateAsAvailability(i, v));
      let index = ui.divText((asAvailabilityStatuses.length - i).toString());
      index.style.marginLeft = '7px';
      asPattern.append(
        ui.divV([
          asAvailabilityStatuses[i],
          index
        ])
      );
    }
    const five = ui.divText("5'");
    five.style.marginTop = '10px';
    asPattern.append(five);
  }

  function updateAsModification() {
    asModificationItems.innerHTML = '';
    asPtoStatuses = (sequenceLength.value > asPtoStatuses.length) ?
      asPtoStatuses.concat(Array(sequenceLength.value - asPtoStatuses.length).fill(fullyPto)) :
      asPtoStatuses.slice(asPtoStatuses.length - sequenceLength.value);
    asBaseStatuses = (sequenceLength.value > asBaseStatuses.length) ?
      asBaseStatuses.concat(Array(sequenceLength.value - asBaseStatuses.length).fill(sequenceBase)) :
      asBaseStatuses.slice(asBaseStatuses.length - sequenceLength.value);
    for (let i = 0; i < asLength; i++) {
      asPtoStatuses[i] = ui.boolInput('', asPtoStatuses[i].value);
      asBaseStatuses[i] = ui.choiceInput('', asBaseStatuses[i].value, baseChoices);
      asModificationItems.append(
        ui.divH([
          ui.block25([ui.label((i + 1).toString())])!,
          ui.block50([asBaseStatuses[i]])!,
          ui.block25([asPtoStatuses[i]])!
        ], {style: {alignItems: "center"}})
      );
    }
  }

  function updateSsModification() {
    ssModificationItems.innerHTML = '';
    ssPtoStatuses = (sequenceLength.value > ssPtoStatuses.length) ?
      ssPtoStatuses.concat(Array(sequenceLength.value - ssPtoStatuses.length).fill(fullyPto)) :
      ssPtoStatuses.slice(ssPtoStatuses.length - sequenceLength.value);
    ssBaseStatuses = (sequenceLength.value > ssBaseStatuses.length) ?
      ssBaseStatuses.concat(Array(sequenceLength.value - ssBaseStatuses.length).fill(sequenceBase)) :
      ssBaseStatuses.slice(ssBaseStatuses.length - sequenceLength.value);
    for (let i = 0; i < ssLength; i++) {
      ssPtoStatuses[i] = ui.boolInput('', ssPtoStatuses[i].value);
      ssBaseStatuses[i] = ui.choiceInput('', ssBaseStatuses[i].value, baseChoices);
      ssModificationItems.append(
        ui.divH([
          ui.block25([ui.label((i + 1).toString())])!,
          ui.block50([ssBaseStatuses[i]])!,
          ui.block25([ssPtoStatuses[i]])!
        ], {style: {alignItems: "center"}})
      );
    }
  }

  function updateUiForNewSequenceLength() {
    if (sequenceLength.value < maximalValidSequenceLength) {
      updateSsModification();
      updateAsModification();
      updateSsPattern();
      updateAsPattern();
    } else {
      ui.dialog('Sequence length should be less than ' + maximalValidSequenceLength.toString() + ' due to UI constrains')
        .add(ui.divText('Please change sequence length in order to define new pattern'))
        .show();
    }
  }

  function updatePtoAllAtOnce(newPtoValue: boolean) {
    for (let i = 0; i < ssPtoStatuses.length; i++) {
      ssPtoStatuses[i].value = newPtoValue;
    }
    for (let i = 0; i < asPtoStatuses.length; i++) {
      asPtoStatuses[i].value = newPtoValue;
    }
  }

  function updateBasisAllAtOnce(newBasisValue: string) {
    for (let i = 0; i < ssBaseStatuses.length; i++) {
      ssBaseStatuses[i].value = newBasisValue;
    }
    for (let i = 0; i < asBaseStatuses.length; i++) {
      asBaseStatuses[i].value = newBasisValue;
    }
  }

  function updateSsAvailability(indexOfClickedCheckbox: number, isClickedCheckboxChecked: boolean) {
    if (isClickedCheckboxChecked) {
      ssLength = sequenceLength.value - indexOfClickedCheckbox;
      for (let i = indexOfClickedCheckbox; i < sequenceLength.value; i++) {
        ssAvailabilityStatuses[i] = ui.boolInput('', true);
      }
    } else {
      ssLength = sequenceLength.value - indexOfClickedCheckbox - 1;
      for (let i = 0; i < indexOfClickedCheckbox; i++) {
        ssAvailabilityStatuses[i] = ui.boolInput('', false);
      }
    }
    updateSsPattern();
    updateSsModification();
  }

  function updateAsAvailability(indexOfClickedCheckbox: number, isClickedCheckboxChecked: boolean) {
    if (isClickedCheckboxChecked) {
      asLength = sequenceLength.value - indexOfClickedCheckbox;
      for (let i = indexOfClickedCheckbox; i < sequenceLength.value; i++) {
        asAvailabilityStatuses[i] = ui.boolInput('', true);
      }
    } else {
      asLength = sequenceLength.value - indexOfClickedCheckbox - 1;
      for (let i = 0; i < indexOfClickedCheckbox; i++) {
        asAvailabilityStatuses[i] = ui.boolInput('', false);
      }
    }
    updateAsPattern();
    updateAsModification();
  }

  function addColumnsFields() {
    chooseSsColumnDiv.innerHTML = '';
    chooseAsColumnDiv.innerHTML = '';
    let chosenInputSsColumn = ui.choiceInput('SS Column', '', grok.shell.table(tables.value).columns.names());
    let chosenInputAsColumn = ui.choiceInput('AS Column', '', grok.shell.table(tables.value).columns.names());
    chooseSsColumnDiv.append(chosenInputSsColumn.root);
    if (createAsStrand.value) chooseAsColumnDiv.append(chosenInputAsColumn.root);
    patternDesignSection.append(
      ui.button('Convert Sequences', () => convertSequence(chosenInputSsColumn.value, chosenInputAsColumn.value))
    );
  }

  function convertSequence(chosenInputSsColumn: string, chosenInputAsColumn: string) {
    let count: number = -1;
    grok.shell.table(tables.value).columns.addNewString('Axolabs ' + chosenInputSsColumn).init((i: number) => {
      count = -1;
      return grok.shell.table(tables.value).columns.byName(chosenInputSsColumn).get(i).replace(/[AUGC]/g, function (x: string) {
        count++;
        let ind = axolabsMap["RNA"]["symbols"].indexOf(x);
        let v = axolabsMap[ssBaseStatuses[count].value]["symbols"][ind];
        return (ssPtoStatuses[count].value) ? v + 'ps' : v;
      })
    });
    if (createAsStrand.value)
      grok.shell.table(tables.value).columns.addNewString('Axolabs ' + chosenInputAsColumn).init((i: number) => {
        count = -1;
        return grok.shell.table(tables.value).columns.byName(chosenInputAsColumn).get(i).replace(/[AUGC]/g, function (x: string) {
          count++;
          let ind = axolabsMap["RNA"]["symbols"].indexOf(x);
          let v = axolabsMap[asBaseStatuses[count].value]["symbols"][ind];
          return (asPtoStatuses[count].value) ? v + 'ps' : v;
        });
      });
  }

  let ssModificationItems = ui.div([]),
    asModificationItems = ui.div([]),
    ssPattern = ui.divH([]),
    asPattern = ui.divH([]),
    chooseAsColumnDiv = ui.div([]),
    chooseSsColumnDiv = ui.div([]);

  let ssAvailabilityStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultAvailability)),
    asAvailabilityStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultAvailability)),
    ssBaseStatuses = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices)),
    asBaseStatuses = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices)),
    ssPtoStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto)),
    asPtoStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto));

  let sequenceLength = ui.intInput('Sequence Length', defaultSequenceLength, () => {
    if (sequenceLength.value > oldSequenceLength) {
      ssLength += sequenceLength.value - oldSequenceLength;
      asLength += sequenceLength.value - oldSequenceLength;
    }
    if (asLength > sequenceLength.value) asLength = sequenceLength.value;
    if (ssLength > sequenceLength.value) ssLength = sequenceLength.value;
    updateUiForNewSequenceLength();
    oldSequenceLength = sequenceLength.value;
  });
  let oldSequenceLength: number = sequenceLength.value,
    ssLength: number = defaultSequenceLength,
    asLength: number = defaultSequenceLength;

  let tables = ui.choiceInput('Tables', '', grok.shell.tableNames, (name: string) => addColumnsFields());

  let sequenceBase = ui.choiceInput('Sequence Basis', defaultBase, baseChoices, (v: string) => updateBasisAllAtOnce(v));

  let fullyPto = ui.boolInput('Fully PTO', defaultPto, (v: boolean) => updatePtoAllAtOnce(v));

  let createAsStrand = ui.boolInput('Create AS Strand', true, (v: boolean) => {
    asModificationSection.hidden = (!v);
    asPattern.hidden = (!v);
    chooseAsColumnDiv.hidden = (!v);
  });

  let applyExistingDesign = ui.choiceInput('Apply Existing Pattern', '', ['Var-3A97'], () => grok.shell.info('Coming soon'));

  let threeModification = ui.stringInput("Addidional 3' Modification", "", (v: string) => grok.shell.info('Coming soon'));
  let fiveModification = ui.stringInput("Addidional 5' Modification", "", (v: string) => grok.shell.info('Coming soon'));

  sequenceLength.root.append(ui.button('Specify Duplex', () => {
    ui.dialog('Axolabs Pattern')
      .add(
        ui.span([
          drawAxolabsPattern(
            applyExistingDesign.value,
            ssBaseStatuses.slice(0, ssLength).map((e) => e.value),
            asBaseStatuses.slice(0, asLength).map((e) => e.value),
            ssPtoStatuses.slice(0, ssLength).map((e) => e.value),
            asPtoStatuses.slice(0, asLength).map((e) => e.value)
          )
        ])
      )
      .show();
  }));
  updateUiForNewSequenceLength();

  let patternDesignSection = ui.divV([
    ui.h1('Pattern Design'),
    ui.divH([
      ui.div([
        tables.root,
        chooseSsColumnDiv,
        chooseAsColumnDiv,
        applyExistingDesign.root,
        sequenceLength.root,
        sequenceBase.root,
        fullyPto.root,
        createAsStrand.root
      ], 'ui-form'),
      ui.inputs([
        threeModification,
        fiveModification
      ], {})
    ], {style: {flexWrap: 'wrap'}}),
    ui.block([
      ssPattern,
      asPattern
    ], {style: {overflowX: 'scroll'}}),
    ui.divH([
      ui.button('Define Pattern', () => grok.shell.info('Coming soon')),
      ui.button('Save Pattern', () => grok.shell.info('Coming soon'))
    ])
  ]);

  let ssModificationSection = ui.divV([
    ui.h1('Sense Strand'),
    ui.divH([
      ui.block25([ui.divText('#')])!,
      ui.block50([ui.divText('Modification')])!,
      ui.block25([ui.divText('PTO')])!
    ]),
    ssModificationItems
  ], {style: {marginLeft: "30px", width: "250px"}});

  let asModificationSection = ui.divV([
    ui.h1('Antisense Strand'),
    ui.divH([
      ui.block25([ui.divText('#')])!,
      ui.block50([ui.divText('Modification')])!,
      ui.block25([ui.divText('PTO')])!
    ]),
    asModificationItems
  ], {style: {marginLeft: "30px", width: "250px"}});

  let appAxolabsDescription = ui.info(
    [
      ui.divText("\n How to define new pattern:",{style:{'font-weight':'bolder'}}),
      ui.divText("1. Choose table and columns with sense and antisense strands"),
      ui.divText("2. Choose lengths of both strands by editing checkboxes below"),
      ui.divText("3. Choose basis and PTO status for each nucleotide"),
      ui.divText("4. Set additional modifications for sequence edges"),
      ui.divText("5. Press 'Convert Sequences' button"),
      ui.divText("This will add the result column(s) to the right of the table"),
    ], 'Create and apply Axolabs translation patterns.'
  );

  let splitter = ui.splitH([
    patternDesignSection,
    ui.divH([
      ssModificationSection,
      asModificationSection
    ])
  ], {style: {height: 'inherit'}});

  $(splitter).children('.ui-box').removeClass('ui-box').addClass('ui-div');
  $(splitter).children('.ui-div:nth-child(2)').css('flex-shrink','0');

  return ui.block([
    appAxolabsDescription,
    ui.block([
      splitter
    ], {style: {height: '100%'}})
  ], {style: {padding: "20px 0", overflowY: "hidden!important"}}); // will not allow to scroll over page
}