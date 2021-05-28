/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/*
SS - sense strand of DNA;
AS - antisense strand of DNA;
base - indicates how to translate input nucleotide;
PTO - indicates whether oligonucleotide is phosphorothioated (ps linkage);
Pattern design - section of dialog for changing length of AS and SS by editing corresponding checkboxes;
SS/AS Modification - sections of dialog for changing base and PTO statuses of one nucleotide at once: for SS or for AS;
*/

const baseChoices: string[] = ["2'OMe RNA (mX)", "2'Fluoro RNA (fX)", "2'OMe U (mU)", "UNA", "GNA", "InvAb", "dX"];
const defaultBase: string = "2'OMe RNA (mX)";
const defaultPto: boolean = false;
const defaultAvailability: boolean = true;
const defaultSequenceLength: number = 23;
const maximalValidSequenceLength: number = 30;

//name: defineAxolabsPattern
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function defineAxolabsPattern(nucleotides: string) {

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
      ssPattern.append(
        ui.divV([
          ui.divText((i + 1).toString()),
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
      asPattern.append(
        ui.divV([
          asAvailabilityStatuses[i],
          ui.divText((asAvailabilityStatuses.length - i).toString())
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

  let ssModificationItems = ui.div([]),
    asModificationItems = ui.div([]),
    ssPattern = ui.divH([]),
    asPattern = ui.divH([]);

  let ssAvailabilityStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultAvailability)),
    asAvailabilityStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultAvailability)),
    ssBaseStatuses = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices)),
    asBaseStatuses = Array(defaultSequenceLength).fill(ui.choiceInput('', defaultBase, baseChoices)),
    ssPtoStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto)),
    asPtoStatuses = Array(defaultSequenceLength).fill(ui.boolInput('', defaultPto));

  let sequenceLength = ui.intInput('Enter sequence length', defaultSequenceLength, () => {
    asLength = (asLength > sequenceLength.value) ? sequenceLength.value : asLength;
    ssLength = (ssLength > sequenceLength.value) ? sequenceLength.value : ssLength;
    updateUiForNewSequenceLength();
  });
  let ssLength: number = defaultSequenceLength,
    asLength: number = defaultSequenceLength;

  let sequenceBase = ui.choiceInput('Sequence basis', defaultBase, baseChoices, (v: string) => updateBasisAllAtOnce(v));

  let fullyPto = ui.boolInput('Fully PTO', defaultPto, (v: boolean) => updatePtoAllAtOnce(v));

  let createAsStrand = ui.boolInput('Create AS strand', true, (v: boolean) => {
    asModificationSection.hidden = (!v);
    asPattern.hidden = (!v);
  });

  let threeModification = ui.stringInput("Addidional 3' modification", "", (v: string) => grok.shell.info('Coming soon'));
  let fiveModification = ui.stringInput("Addidional 5' modification", "", (v: string) => grok.shell.info('Coming soon'));

  updateUiForNewSequenceLength();

  let patternDesignSection = ui.divV([
    ui.h1('Pattern Design'),
    ui.divH([
      sequenceLength.root,
      ui.button('Specify Duplex', () => grok.shell.info('Coming soon'))
    ]),
    sequenceBase.root,
    fullyPto.root,
    createAsStrand.root,
    ssPattern,
    asPattern,
    threeModification.root,
    fiveModification.root,
    ui.divH([
      ui.button('Define Pattern', () => grok.shell.info('Coming soon')),
      ui.button('Convert Sequences', () => grok.shell.info('Coming soon')),
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

  ui.dialog('Axolabs')
    .add(
      ui.divH([
        patternDesignSection,
        ssModificationSection,
        asModificationSection
      ])
    )
    .showModal(true);
}