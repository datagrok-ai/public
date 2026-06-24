import {before, category, test, expect} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {BoltzService} from '../utils/boltz-api-service';
import {BOLTZ_API_KEY_PARAM} from '../utils/boltz-api-constants';
import {
  structureAndBindingRequest, admeRequest,
  smallMoleculeDesignRequest, smallMoleculeLibraryScreenRequest,
  proteinDesignRequest, proteinLibraryScreenRequest,
} from './fixtures';

const TIMEOUT = 60 * 60 * 1000;

category('Hosted API', () => {
  let boltz: BoltzService;

  before(async () => {
    boltz = BoltzService.getInstance();
    const apiKey = (await _package.getCredentials())?.parameters?.[BOLTZ_API_KEY_PARAM];
    await boltz.init(apiKey);
  });

  test('structureAndBinding', async () => {
    const result = await boltz.predictStructureAndBinding(structureAndBindingRequest);
    expect(result.status, 'succeeded');
    expect(result.output!.all_sample_results.length >= 1, true);
    expect(typeof result.output!.all_sample_results[0].metrics.structure_confidence, 'number');
    expect(typeof result.output!.binding_metrics!.binding_confidence, 'number');
  }, {timeout: TIMEOUT});

  test('adme', async () => {
    const result = await boltz.predictAdme(admeRequest);
    expect(result.status, 'succeeded');
    expect(result.output!.molecules.length, admeRequest.input.molecules.length);
    const aspirin = result.output!.molecules[0];
    expect(aspirin.external_id, 'aspirin');
    expect(aspirin.status, 'succeeded');
  }, {timeout: TIMEOUT});

  test('smallMolecule.design', async () => {
    const results = await boltz.designSmallMolecules(smallMoleculeDesignRequest);
    expect(results.length > 0, true);
    expect(typeof results[0].smiles, 'string');
    expect(typeof results[0].metrics.binding_confidence, 'number');
  }, {timeout: TIMEOUT});

  test('smallMolecule.libraryScreen', async () => {
    const results = await boltz.screenSmallMoleculeLibrary(smallMoleculeLibraryScreenRequest);
    expect(results.length > 0, true);
    expect(typeof results[0].smiles, 'string');
    expect(typeof results[0].metrics.binding_confidence, 'number');
  }, {timeout: TIMEOUT});

  test('protein.design', async () => {
    const results = await boltz.designProteins(proteinDesignRequest);
    expect(results.length > 0, true);
    expect(results[0].entities.length > 0, true);
    expect(typeof results[0].metrics.binding_confidence, 'number');
  }, {timeout: TIMEOUT});

  test('protein.libraryScreen', async () => {
    const results = await boltz.screenProteinLibrary(proteinLibraryScreenRequest);
    expect(results.length > 0, true);
    expect(results[0].entities.length > 0, true);
    expect(typeof results[0].metrics.binding_confidence, 'number');
  }, {timeout: TIMEOUT});
});
