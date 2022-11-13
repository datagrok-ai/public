/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {FastaFileHandler} from '@datagrok-libraries/bio/src/utils/fasta-handler';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

category('monomer lib', () => { 
  //   test('monomerManager', async() => {
  //   const df: DG.DataFrame = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/test.csv'));
  //   grok.shell.addTableView(df);
  //   await grok.functions.call('Helm:monomerManager', {value: 'HELMCoreLibrary.json'});
  //   const checkName = '2-Chloroadenine';
  //   let flag = false;
  //   const types = Object.keys(org.helm.webeditor.monomerTypeList());
  //   const monomers: any = [];
  //   for (var i = 0; i < types.length; i++) {
  //     //@ts-ignore
  //     monomers.push(new scil.helm.Monomers.getMonomerSet(types[i]));
  //     Object.keys(monomers[i]).forEach(k => {
  //       if (monomers[i][k].n == checkName){
  //         flag = true;
  //       }
  //     });
  //   }
  //   expect(flag, true);
  // });
});
