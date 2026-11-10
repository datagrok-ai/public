/**
 * A stack overflow out of a Crux call counts as a trap (crux-api.ts `isWasmTrap`): V8 unwinds a wasm stack overflow as
 * `RangeError: Maximum call stack size exceeded` without restoring the instance's stack pointer, so after a few the
 * instance reads and writes out of bounds. The worker jobs rethrow it, the pool replaces the workers, and on the main
 * thread Crux is off for the session, as after a trap. (Crux's recursive writers overflow on molecules of a few
 * thousand atoms: a 20,000-carbon chain.)
 */
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine, cruxReady, isCruxRequested} from '../engine/chem-engine';
import {loadCrux, noteCruxError, resetCruxLoaderForTests} from '../engine/crux-module';
import {CruxMacromoleculeAssembler} from '../engine/biopolymer';
import {isWasmTrap, throwIfTrap} from '../crux/crux-api';
import {getCruxService} from '../crux/crux-searches';
import {inEngine} from './engine-utils';

const OVERFLOW = 'Maximum call stack size exceeded';

category('chemical engine: traps', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('a stack overflow is a trap: the worker jobs rethrow it, the pool restarts', async () => {
    expect(isWasmTrap(new RangeError(OVERFLOW)), true, 'a stack overflow');
    expect(isWasmTrap(Object.assign(new Error('too much recursion'), {name: 'InternalError'})), true,
      'Firefox\'s stack overflow');
    expect(isWasmTrap(new WebAssembly.RuntimeError('unreachable')), true, 'a trap');
    expect(isWasmTrap(new RangeError('Invalid array length')), false, 'another RangeError');
    expect(isWasmTrap(new Error(OVERFLOW)), false, 'not a RangeError');
    let rethrown = false;
    try {
      throwIfTrap(new RangeError(OVERFLOW));
    } catch {
      rethrown = true;
    }
    expect(rethrown, true, 'throwIfTrap rethrows a stack overflow (a job does not turn it into "no result")');

    await inEngine(ChemicalEngine.Crux, async () => {
      const service = await getCruxService();
      // as the worker reports one (crux.worker.ts posts the RangeError's message)
      let generation = service.generation;
      await service.failJobForTests(OVERFLOW).catch(() => {});
      expect(service.generation, generation + 1, 'a worker\'s stack overflow restarts the workers');
      // a real one where Crux's SMILES writer still recurses (a no-op test once it does not)
      generation = service.generation;
      let failure = '';
      try {
        await service.jobInChunks<string[]>('chemCanonicalSmiles', ['C'.repeat(20000)], 50);
      } catch (e: any) {
        failure = String(e?.message ?? e);
      }
      if (failure.includes(OVERFLOW))
        expect(service.generation, generation + 1, 'the real stack overflow restarted the workers');
      expect(isCruxRequested(), true, 'a worker\'s failure leaves Crux mode on');
    });
  }, {timeout: 60000});

  test('a stack overflow on the main thread leaves RDKit for the session', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      try {
        noteCruxError(new RangeError('Invalid array length'));
        expect(isCruxRequested(), true, 'another RangeError leaves Crux on');
        noteCruxError(new RangeError(OVERFLOW));
        expect(isCruxRequested(), false, 'after a stack overflow: RDKit');
      } finally {
        resetCruxLoaderForTests(null); // the next session's Crux
      }
      expect(await cruxReady(), true, 'Crux back after the reset');
    });
  });

  test('after a trap on the main thread the assembler\'s synchronous calls refuse', async () => {
    // as its async calls do (they load the module, which fails for the session): a trapped instance may be broken
    await inEngine(ChemicalEngine.Crux, async () => {
      const crux = await loadCrux('test');
      const mol = crux.get_mol('CCO')!;
      const molblock = mol.get_v3Kmolblock();
      mol.delete();
      const assembler = new CruxMacromoleculeAssembler();
      const v2000 = await assembler.toV2000(molblock); // loads the main thread's module for the synchronous twins
      expect(assembler.toV2000Sync(molblock), v2000, 'the synchronous twin before the trap');
      try {
        noteCruxError(new WebAssembly.RuntimeError('simulated: unreachable'));
        let refused = '';
        try {
          assembler.toV2000Sync(molblock);
        } catch (e: any) {
          refused = String(e?.message ?? e);
        }
        expect(refused.includes('trapped'), true, `toV2000Sync after the trap: ${refused || 'answered'}`);
        let asyncRefused = false;
        await assembler.toV2000(molblock).catch(() => asyncRefused = true);
        expect(asyncRefused, true, 'toV2000 after the trap');
      } finally {
        resetCruxLoaderForTests(null); // the next session's Crux
      }
      expect(await cruxReady(), true, 'Crux back after the reset');
    });
  });
});
