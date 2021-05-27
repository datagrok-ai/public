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
const defaultSequenceLength: number = 23;

//name: defineAxolabsPattern
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function defineAxolabsPattern(nucleotides: string) {

  function updateSsPattern(sequenceLength: number) {
    ssPattern.innerHTML = '';
    const ss5 = ui.divText("SS-5'");
    ss5.style.marginTop = '26px';
    ssPattern.append(ss5);
    for (let i = 0; i < sequenceLength; i++) {
      ssPattern.append(
        ui.divV([
          ui.divText((i + 1).toString()),
          ui.boolInput('', true, (v: boolean) => grok.shell.info(v)).root
        ])
      );
    }
    const three = ui.divText("3'");
    three.style.marginTop = '26px';
    ssPattern.append(three);
  }

  function updateAsPattern(sequenceLength: number) {
    asPattern.innerHTML = '';
    const as3 = ui.divText("AS-3'");
    as3.style.marginTop = '10px';
    asPattern.append(as3);
    for (let i = 0; i < sequenceLength; i++) {
      asPattern.append(
        ui.divV([
          ui.boolInput('', true, (v: boolean) => grok.shell.info(v)).root,
          ui.divText((i + 1).toString())
        ])
      );
    }
    const five = ui.divText("5'");
    five.style.marginTop = '10px';
    asPattern.append(five);
  }

  function updateModification(section: HTMLElement, sequenceLength: number, baseStatuses: string[], ptoStatuses: boolean[]) {
    section.innerHTML = '';

    for (let i = 0; i < sequenceLength; i++) {
      section.append(
        ui.divH([
          ui.block25([ui.label((i + 1).toString())])!,
          ui.block50([ui.choiceInput('', baseStatuses[i], baseChoices, (v: string) => {
            baseStatuses[i] = v
          }).root])!,
          ui.block25([
            ui.div([
              ui.boolInput('', ptoStatuses[i], (v: boolean) => {
                $('.modif .ui-input-bool .ui-input-editor').on('click', function() {
                  let ind = parseInt($(this).parent().parent().attr('class')!.slice(7));
                  ptoStatuses[ind] = v
                });
              }).root
            ], 'myinfo-' + (i + 1).toString())
          ])!
        ], {style: {alignItems: "center"}})
      )
    }
  }

  function updateUiForNewSequenceLength(sequenceLength: number) {
    updateModification(asModification, sequenceLength, asBaseStatuses, asPtoStatuses);
    updateModification(ssModification, sequenceLength, ssBaseStatuses, ssPtoStatuses);
    updateSsPattern(sequenceLength);
    updateAsPattern(sequenceLength);
  }

  let ssModification = ui.div([]),
    asModification = ui.div(),
    ssPattern = ui.divH([]),
    asPattern = ui.divH([]);

  let ssBaseStatuses: string[] = Array(defaultSequenceLength).fill(defaultBase),
    asBaseStatuses: string[] = Array(defaultSequenceLength).fill(defaultBase),
    ssPtoStatuses: boolean[] = Array(defaultSequenceLength).fill(defaultPto),
    asPtoStatuses: boolean[] = Array(defaultSequenceLength).fill(defaultPto);

  let sequenceLength = ui.intInput('Enter sequence length', defaultSequenceLength, (v: number) => updateUiForNewSequenceLength(v));

  let sequenceBase = ui.choiceInput('Sequence basis', defaultBase, baseChoices, (v: string) => {
    updateModification(asPattern, sequenceLength.value, asBaseStatuses, asPtoStatuses);
    updateModification(ssPattern, sequenceLength.value, ssBaseStatuses, ssPtoStatuses);
  });

  let fullyPto = ui.boolInput('Fully PTO', defaultPto, (v: boolean) => {grok.shell.info('Coming soon')});

  let createAsStrand = ui.boolInput('Create AS strand', true, (v: boolean) => {
    asModificationSection.hidden = (!v);
    asPattern.hidden = (!v);
  });

  let threeModification = ui.stringInput("Addidional 3' modification", "", (v: string) => grok.shell.info('Coming soon'));
  let fiveModification = ui.stringInput("Addidional 5' modification", "", (v: string) => grok.shell.info('Coming soon'));

  updateUiForNewSequenceLength(defaultSequenceLength);

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
    ssModification,
    threeModification.root,
    fiveModification.root
  ], {style: {marginLeft: "30px", width: "250px"}});

  let asModificationSection = ui.divV([
    ui.h1('Antisense Strand'),
    ui.divH([
      ui.block25([ui.divText('#')])!,
      ui.block50([ui.divText('Modification')])!,
      ui.block25([ui.divText('PTO')])!
    ]),
    asModification
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

