/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

//name: OligoBatchCalculator
//tags: app
export function OligoBatchCalculator() {

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let text1 = ui.divText('Enter Yield Amount & Units');
  let text2 = ui.divText('Search Modifications');

  let enter1 = ui.stringInput('', '');
  let enter2 = ui.stringInput('', '');
  let choose = ui.choiceInput('', 'Select', []);
  let analyzeAsSingleStrand = ui.button('Analyze As Single Strand', () => grok.shell.info('Coming soon'));
  let analyzeAsDuplexStrand = ui.button('Analyze As Duplex Strand', () => grok.shell.info('Coming soon'));
  let clearSequences = ui.button('Clear Sequences', () => grok.shell.info('Coming soon'));
  let addModifications = ui.button('Add Modifications', () => grok.shell.info('Coming soon'));
  let threeMod = ui.boolInput("3' MOD", false);
  let internal = ui.boolInput('INTERNAL', false);
  let fiveMod = ui.boolInput("5' MOD", false);
  let inputSequenceField = ui.textInput("", "", async (seq: string) => {});

  let tableRows = [];
  for (let i = 0; i < 10; i++) {
    // @ts-ignore
    tableRows.push({'key': i, 'value': i*19})
  }

  grok.shell.newView('Sequence Translator', [
    ui.divH([
      ui.divV([
        text1,
        ui.divH([
          enter1.root,
          choose.root
        ])
      ]),
      analyzeAsSingleStrand,
      analyzeAsDuplexStrand,
      clearSequences,
      addModifications,
      ui.divV([
        text2,
        enter2.root
      ]),
      threeMod.root,
      internal.root,
      fiveMod.root
    ]),
    ui.block([
      ui.div([
        ui.h1('Input sequence'),
        ui.div([
          inputSequenceField.root
        ],'input-base')
      ], 'sequenceInput'),

    ]),
    ui.button('Load Data To Excel', () => {}),
    DG.HtmlTable.create(
      tableRows,
      (item: {key: string; value: string;}) => [item.key, item.value],
      ['Item#', 'Sequence', 'Length', 'OD 250', 'nmole', 'Mass', 'nmole/OD', 'μg/OD', 'MW', 'G%', 'GC%', 'Ext. Coefficient']
    ).root
  ]);

  $('.sequence')
    .children().css('padding','5px 0');
  $('.sequenceInput .input-base').css('margin','0');
  $('.sequenceInput textarea')
    .css('resize','none')
    .css('min-height','50px')
    .css('width','100%');
  $('.sequenceInput select')
    .css('width','100%');
}
