import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {awaitCheck, category, test, testViewer} from './_timed-test';

const _tsLog = (msg: string): void => console.log(`[${new Date().toISOString()}] ${msg}`);

async function awaitChemSearchViewer(viewer: DG.Viewer): Promise<void> {
  _tsLog(`[VIEWER-TEST] awaitChemSearchViewer: entering for type=${viewer.type}`);
  await awaitCheck(()=> {
    const className = viewer.type === 'Chem Similarity Search' ? 'chem-viewer-grid' : 'chem-diversity-search';
    const el = document.getElementsByClassName(className);
    const errEl = document.getElementsByClassName('chem-malformed-molecule-error');
    return (el.length && el[0].children.length === 12) || !!errEl.length;
  }, 'Viewer hasn\'t been rendered', 3000);
  _tsLog(`[VIEWER-TEST] awaitChemSearchViewer: done for type=${viewer.type}`);
}

async function awaitScaffoldTree(viewer: DG.Viewer): Promise<void> {
  await awaitCheck(() => document.getElementsByClassName('mol-host').length > 0,
    'scaffold tree has not been generated', 150000);
}

const awaitViewers: {[key: string]: (viewer: DG.Viewer) => Promise<void>} = {
  'Chem Similarity Search': awaitChemSearchViewer,
  'Chem Diversity Search': awaitChemSearchViewer,
  'Scaffold Tree': awaitScaffoldTree,
}


category('viewers', () => {
  const df = grok.data.demo.molecules(15);
  const viewers = DG.Func.find({ package: 'Chem', meta: {role: DG.FUNC_TYPES.VIEWER} }).reduce<string[]>((result, f) => {
    const name = f.friendlyName;
    if (name !== 'Scaffold Tree')
      result.push(name);
    return result;
  }, []);

  _tsLog(`[VIEWER-TEST] viewers category: discovered ${viewers.length} viewers: ${JSON.stringify(viewers)}`);

  for (const v of viewers) {
    test(v, async () => {
      const t0 = performance.now();
      _tsLog(`[VIEWER-TEST] "${v}": entering test, calling testViewer with df.rows=${df.rowCount}, ` +
        `awaitViewer=${!!awaitViewers[v]}`);
      try {
        //@ts-ignore
        await testViewer(v, df.clone(), {
          detectSemanticTypes: true,
          readOnly: false,
          arbitraryDfTest: false,
          awaitViewer: awaitViewers[v]
        });
        _tsLog(`[VIEWER-TEST] "${v}": testViewer returned after ${(performance.now() - t0).toFixed(0)} ms`);
      } catch (e: any) {
        const msg = e instanceof Error ? e.message : String(e);
        _tsLog(`[VIEWER-TEST] "${v}": testViewer threw after ${(performance.now() - t0).toFixed(0)} ms: ${msg}`);
        throw e;
      }
    }, {timeout: 90000});
  }
});