import { after, awaitCheck, before, category, delay, expect, test } from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { _testSetMolfile, _testSetSmarts, _testSetSmiles } from './ketcher-utils';
import { KetcherSketcher } from '../ketcher';

/** A Ketcher sketcher in a dialog of its own, once its Ketcher has started. */
async function openKetcher(): Promise<{ketcher: any, dialog: DG.Dialog}> {
  const func = DG.Func.find({meta: {role: DG.FUNC_TYPES.MOLECULE_SKETCHER}, name: 'ketcherSketcher'})[0];
  grok.chem.currentSketcherType = func.friendlyName;
  const s = new grok.chem.Sketcher();
  const dialog = ui.dialog().add(s).show();
  await s.sketcherReady();
  await awaitCheck(() => (s.sketcher as KetcherSketcher)?._sketcher !== null, 'Ketcher did not start', 20000);
  return {ketcher: (s.sketcher as KetcherSketcher)._sketcher!, dialog};
}


category('ketcher', async () => {
  let previousSketcherType: string | undefined;

  before(async () => {
    previousSketcherType = grok.chem.currentSketcherType;
  });

  after(async () => {
    if (previousSketcherType !== undefined)
      grok.chem.currentSketcherType = previousSketcherType;
  });

  test('setSmiles', async () => {
    await _testSetSmiles();
  });

  test('setMolfile', async () => {
    await _testSetMolfile();
  });

  test('setSmarts', async () => {
    await _testSetSmarts();
  });

  // Indigo's struct service leaves a conversion unanswered when another meets it; the page's queue of conversions
  // stalled for good behind one before (crux-sketch spike query-roundtrip)
  test('a conversion Indigo never answers holds up no later one', async () => {
    const inTurn = (KetcherSketcher as any)._inTurn as <T>(c: () => Promise<T>) => Promise<T>;
    const unanswered = inTurn(() => new Promise<string>(() => {}));
    const next = inTurn(() => Promise.resolve('answered'));
    let given = '';
    await unanswered.catch((e: any) => given = String(e?.message ?? e));
    expect(given, 'Indigo did not answer', 'the unanswered conversion');
    expect(await next, 'answered', 'the conversion after it');
  }, {timeout: 40000});

  // why one went unanswered: ketcher-standalone's one Indigo worker of the page answers through the IndigoService
  // created last, and a reply drops another request of its kind in flight (crux-sketch spike query-roundtrip, K3)
  const everyRequest = 'Indigo answers every request of the page\'s Ketchers: two at once, and one of a Ketcher ' +
    'mounted before';
  test(everyRequest, async () => {
    const smiles = 'chemical/x-daylight-smiles' as any;
    const answer = (ketcher: any, input: string): Promise<string> => Promise.race([
      ketcher.indigo.convert(input, {outputFormat: smiles}).then((r: any) => String(r.struct).trim(),
        (e: any) => `failed: ${e?.message ?? e}`),
      delay(10000).then(() => 'unanswered')]);
    const first = await openKetcher();
    const second = {dialog: null as DG.Dialog | null};
    try {
      const [ethanol, ethylamine] = await Promise.all([answer(first.ketcher, 'CCO'), answer(first.ketcher, 'CCN')]);
      expect(ethanol, 'CCO', 'the first of two requests at once');
      expect(ethylamine, 'CCN', 'the second of two requests at once');
      const later = await openKetcher();
      second.dialog = later.dialog;
      expect(await answer(first.ketcher, 'CCS'), 'CCS', 'a request of the Ketcher mounted before');
      expect(await answer(later.ketcher, 'CCC'), 'CCC', 'a request of the Ketcher mounted last');
    } finally {
      second.dialog?.close();
      first.dialog.close();
    }
  }, {timeout: 90000});

  // ketcher-core's Ketcher never ended its subscription to the page's one SettingsService, so every Ketcher
  // mounted stayed reachable, and a page that had mounted some 50 to 150 filled its heap (crux-sketch spike
  // query-roundtrip, K4)
  test('a Ketcher closed or paused leaves no subscription on the page\'s settings service', async () => {
    // the settings service the Ketchers share (ketcher-core's singleton), as each Ketcher holds it
    const services = new Set<any>();
    const opened = async (): Promise<{ketcher: any, dialog: DG.Dialog}> => {
      const k = await openKetcher();
      services.add(k.ketcher.settingsService);
      return k;
    };
    const listeners = (): number => Array.from(services)
      .reduce((n, s) => n + (s?.emitter?.listenerCount?.('settings:changed') ?? 0), 0);
    const closed = async (n: number): Promise<void> => {
      for (let i = 0; i < n; i++) {
        const {dialog} = await opened();
        dialog.close();
        await delay(300);
      }
    };
    await closed(1);
    const settled = listeners();
    await closed(3);
    expect(listeners(), settled,
      `the settings listeners after three more Ketchers opened and closed (${services.size} settings services)`);
    const paused = await opened();
    const last = await opened();
    last.dialog.close();
    paused.dialog.close();
    await delay(300);
    expect(listeners(), settled, 'the settings listeners after a Ketcher paused by another, and both closed');
  }, {timeout: 120000});

});
