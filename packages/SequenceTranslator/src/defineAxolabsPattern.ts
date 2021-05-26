/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


//name: defineAxolabsPattern
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function defineAxolabsPattern(nucleotides: string) {

  function updateUiForNewSequenceLength(sequenceLength: number) {

    senseStrandModification.innerHTML = '';
    antiSenseStrandModification.innerHTML = '';
    senseStrandPattern.innerHTML = '';
    antiSenseStrandPattern.innerHTML = '';

    senseStrandPattern.append(ss5);
    antiSenseStrandPattern.append(as3);
    for (let i = 0; i < sequenceLength; i++) {
      senseStrandModification.append(
        ui.divH([
          ui.block25([ui.label((i + 1).toString())])!,
          ui.block50([ui.choiceInput('', sequenceBasis.value, choices).root])!,
          ui.block25([ui.boolInput('', false).root])!
        ], {style:{alignItems:"center"}})
      );
      antiSenseStrandModification.append(
        ui.divH([
          ui.block25([ui.label((i + 1).toString())])!,
          ui.block50([ui.choiceInput('', sequenceBasis.value, choices).root])!,
          ui.block25([ui.boolInput('', false).root])!
        ], {style:{alignItems:"center"}})
      );
      senseStrandPattern.append(
        ui.divV([
          ui.divText((i + 1).toString()),
          ui.boolInput('', true).root
        ])
      );
      senseStrandPattern.append(three);
      antiSenseStrandPattern.append(
        ui.divV([
          ui.boolInput('', true).root,
          ui.divText((sequenceLength - i).toString())
        ])
      );
      antiSenseStrandPattern.append(five);
    }
  }

  function updateSequenceBasis(sequenceBasis: string) {
    for (let i = 0; i < 5; i++) {

    }
  }

  const choices = ["2'OMe RNA", "2'Fluoro RNA", "UNA"];
  const defaultSequenceLength = 23;
  const defaultSequenceBasis = "2'OMe RNA";

  const ss5 = ui.divText("SS-5'");
  const as3 = ui.divText("AS-3'");
  const three = ui.divText("3'");
  const five = ui.divText("5'");
  ss5.style.marginTop = '26px';
  as3.style.marginTop = '10px';
  three.style.marginTop = '26px';
  five.style.marginTop = '10px';

  let senseStrandModification = ui.div(),
    antiSenseStrandModification = ui.div(),
    senseStrandPattern = ui.divH([]),
    antiSenseStrandPattern = ui.divH([]);

  let sequenceLength = ui.intInput('Enter sequence length', defaultSequenceLength, (v: number) => updateUiForNewSequenceLength(v));
  let sequenceBasis = ui.choiceInput('Sequence basis', defaultSequenceBasis, choices, (v: string) => updateSequenceBasis(v));
  let fullyPto = ui.boolInput('Fully PTO', false);
  let createAsStrand = ui.boolInput('Create AS strand', true, (v: boolean) => {
    antiSenseStrandSection.hidden = (!v);
    antiSenseStrandPattern.hidden = (!v);
  });
  let threeModification = ui.stringInput("Addidional 3' modification", "", () => grok.shell.info('Coming soon'));
  let fiveModification = ui.stringInput("Addidional 5' modification", "", () => grok.shell.info('Coming soon'));

  updateUiForNewSequenceLength(defaultSequenceLength);

  let defineIndividualDesign = ui.divV([
    ui.h1('Define individual design'),
    ui.divH([
      sequenceLength.root,
      ui.button('Specify duplex', () => grok.shell.info('Coming soon'))
    ]),
    sequenceBasis.root,
    fullyPto.root,
    createAsStrand.root,
    senseStrandPattern,
    antiSenseStrandPattern,
    ui.divH([
      ui.button('Define pattern', () => grok.shell.info('Coming soon')),
      ui.button('Convert sequences', () => grok.shell.info('Coming soon')),
      ui.button('Save pattern', () => grok.shell.info('Coming soon'))
    ])
  ]);

  let senseStrandSection = ui.divV([
    ui.h1('Sense Strand'),
    ui.divH([
      ui.block25([ui.divText('#')])!,
      ui.block50([ui.divText('Modification')])!,
      ui.block25([ui.divText('PTO')])!
    ]),
    senseStrandModification,
    threeModification.root,
    fiveModification.root
  ], {style:{marginLeft: "30px", width: "250px"}});

  let antiSenseStrandSection = ui.divV([
    ui.h1('Antisense Strand'),
    ui.divH([
      ui.block25([ui.divText('#')])!,
      ui.block50([ui.divText('Modification')])!,
      ui.block25([ui.divText('PTO')])!
    ]),
    antiSenseStrandModification
  ], {style:{marginLeft: "30px", width: "250px"}});

  let inp = ui.divH([
    defineIndividualDesign,
    senseStrandSection,
    antiSenseStrandSection
  ]);

  ui.dialog('Axolabs')
    .add(inp)
    .showModal(true);
}