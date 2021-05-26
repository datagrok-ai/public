/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


//name: defineAxolabsPattern
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function defineAxolabsPattern(nucleotides: string) {

  function updateUiForNewSequenceLength(sequenceLength: number) {

    senseStrand.innerHTML = '';
    antiSenseStrand.innerHTML = '';
    defineIndividualDesign.innerHTML = '';

    for (let i = 0; i < sequenceLength; i++) {
      senseStrand.append(
        ui.divH([
          ui.choiceInput((i + 1).toString(), sequenceBasis.value, choices).root,
          ui.boolInput('', true).root
        ])
      );
      antiSenseStrand.append(
        ui.divH([
          ui.choiceInput((i + 1).toString(), sequenceBasis.value, choices).root,
          ui.boolInput('', true).root
        ])
      );
      defineIndividualDesign.append(
        ui.divV([
          ui.divText(i.toString()),
          ui.boolInput('', true).root,
          ui.boolInput('', true).root,
          ui.divText((sequenceLength - i).toString())
        ])
      );
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
  as3.style.marginTop = '20px';
  three.style.marginTop = '26px';
  five.style.marginTop = '20px';

  let senseStrand = ui.div(),
    antiSenseStrand = ui.div(),
    defineIndividualDesign = ui.divH([]);

  let sequenceLength = ui.intInput('Enter sequence length', defaultSequenceLength, (v: number) => updateUiForNewSequenceLength(v));
  let sequenceBasis = ui.choiceInput('Sequence basis', defaultSequenceBasis, choices, (v: string) => updateSequenceBasis(v));
  let fullyPto = ui.boolInput('Fully PTO', false);
  let createAsStrand = ui.boolInput('Create AS strand', true);
  let threeModification = ui.stringInput("Addidional 3' modification", "", () => grok.shell.info('Coming soon'));
  let fiveModification = ui.stringInput("Addidional 3' modification", "", () => grok.shell.info('Coming soon'));

  updateUiForNewSequenceLength(defaultSequenceLength);

  let inp = ui.divH([
    ui.divV([
      ui.h1('Define individual design'),
      ui.divH([
        sequenceLength.root,
        ui.button('Specify duplex', () => grok.shell.info('Coming soon'))
      ]),
      sequenceBasis.root,
      fullyPto.root,
      createAsStrand.root,
      ui.divH([
        ui.divV([ss5, as3]),
        defineIndividualDesign,
        ui.divV([three, five])
      ]),
      ui.divH([
        ui.button('Define pattern', () => grok.shell.info('Coming soon')),
        ui.button('Convert sequences', () => grok.shell.info('Coming soon')),
        ui.button('Save pattern', () => grok.shell.info('Coming soon'))
      ])
    ]),
    ui.divV([
      ui.h1('Sense Strand'),
      ui.divH([ui.divText('#'), ui.divText('Modification'), ui.divText('PTO')]),
      senseStrand,
      threeModification.root,
      fiveModification.root
    ]),
    ui.divV([
      ui.h1('Antisense Strand'),
      ui.divH([ui.divText('#'), ui.divText('Modification'), ui.divText('PTO')]),
      antiSenseStrand
    ])
  ]);

  ui.dialog('Axolabs')
    .add(inp)
    .showModal(true);
}