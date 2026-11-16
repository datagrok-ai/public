/**
 * A stack overflow out of a Crux call counts as a trap (crux-api.ts `isWasmTrap`): V8 unwinds a wasm stack overflow as
 * `RangeError: Maximum call stack size exceeded` without restoring the instance's stack pointer, so after a few the
 * instance reads and writes out of bounds. A worker job that reports one has the pool replace the workers; on the main
 * thread the module recovers inside the call, as after any trap (crux/crux-recover.ts), and Crux stays on. (Crux's
 * recursive writers overflow on molecules of a few thousand atoms: a 20,000-carbon chain.)
 */
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine, isCruxRequested} from '../engine/chem-engine';
import {armCruxTrapForTests, loadCrux, noteCruxError} from '../engine/crux-module';
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

  // a trap never switches Crux off: the main thread's module recovers
  // inside the call (crux/crux-recover.ts), and a call that overflows on a fresh instance too fails alone
  test('a stack overflow on the main thread fails that call alone: Crux stays on', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const crux = await loadCrux('test');
      noteCruxError(new RangeError('Invalid array length'));
      expect(isCruxRequested(), true, 'another RangeError leaves Crux on');
      // a real one where Crux's SMILES writer still recurses (a no-op test once it does not)
      const chain = crux.get_mol('C'.repeat(20000));
      let failure = '';
      try {
        chain?.get_smiles();
      } catch (e: any) {
        failure = String(e?.message ?? e);
        noteCruxError(e);
      } finally {
        chain?.delete();
      }
      expect(failure.includes(OVERFLOW), false, `the overflow reached the caller: ${failure}`);
      expect(isCruxRequested(), true, 'Crux stays on');
      const mol = crux.get_mol('CCO')!;
      try {
        expect(mol.get_smiles(), 'CCO', 'the next call answers');
      } finally {
        mol.delete();
      }
    });
  }, {timeout: 60000});

  test('a trap on the main thread in the assembler\'s synchronous call is recovered from: it answers', async () => {
    // crux-core's real test trap (the next engine call panics), not a stub
    await inEngine(ChemicalEngine.Crux, async () => {
      const crux = await loadCrux('test');
      const mol = crux.get_mol('CCO')!;
      const molblock = mol.get_v3Kmolblock();
      mol.delete();
      const assembler = new CruxMacromoleculeAssembler();
      const v2000 = await assembler.toV2000(molblock); // loads the main thread's module for the synchronous twins
      expect(assembler.toV2000Sync(molblock), v2000, 'the synchronous twin before the trap');
      expect(armCruxTrapForTests(1), true, 'the main thread\'s module is loaded');
      expect(assembler.toV2000Sync(molblock), v2000, 'toV2000Sync, a trap in it');
      expect(isCruxRequested(), true, 'Crux stays on');
      expect(await assembler.toV2000(molblock), v2000, 'toV2000 after the trap');
    });
  });
});
