import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, test, expect, delay, expectFloat, before} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {runStructuralAlertsDetection} from '../panels/structural-alerts';
import {RDMol} from '../rdkit-api';

category('Screening tools benchmarks', () => {
  before(async () => {
    chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
    chemCommonRdKit.initRdKitModuleLocal();
  });

  test('Structural Alerts Benchmark', async () => {
    const alertsDf = DG.DataFrame.fromCsv(await _package.files.readAsText('alert-collection.csv'));
    const ruleSetCol = alertsDf.getCol('rule_set_name');
    const smartsCol = alertsDf.getCol('smarts');
    const ruleIdCol = alertsDf.getCol('rule_id');
    const rdkitModule = chemCommonRdKit.getRdKitModule();
  
    const smartsMap = new Map<string, RDMol>();
    for (let i = 0; i < alertsDf.rowCount; i++)
      smartsMap.set(ruleIdCol.get(i), rdkitModule.get_qmol(smartsCol.get(i)));

    const sarSmall = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/smi10K.csv'));
    const smilesCol = sarSmall.getCol('smiles');
    const ruleSetList = ['BMS', 'Dundee', 'Glaxo', 'Inpharmatica', 'LINT', 'MLSMR', 'PAINS', 'SureChEMBL'];

    DG.time('Structural Alerts', () => {
      runStructuralAlertsDetection(sarSmall, ruleSetList, smilesCol, ruleSetCol, ruleIdCol, smartsMap, rdkitModule);
    });
  });
});