/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function showAppInfo() {
  const textDiv = ui.divText('hello');
  ui.dialog('Dialog name').add(textDiv).show();
}

const appMainDescription = ui.info([
  ui.divText('How to convert one sequence:', {style: {'font-weight': 'bolder'}}),
  ui.divText('Paste sequence into the text field below'),
  ui.divText('\n How to convert many sequences:', {style: {'font-weight': 'bolder'}}),
  ui.divText('1. Drag & drop an Excel or CSV file with sequences into Datagrok'),
  ui.divText('2. Right-click on the column header, then see the \'Convert\' menu'),
  ui.divText('This will add the result column to the right of the table'),
], 'Convert oligonucleotide sequences between Nucleotides, BioSpring, Axolabs, Mermade 12 and GCRS representations.');

