import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {awaitCheck, category, test, testViewer} from '@datagrok-libraries/utils/src/test';

async function awaitChemSearchViewer(viewer: DG.Viewer): Promise<void> {
  await awaitCheck(()=> {
    const className = viewer.type === 'Chem Similarity Search' ? 'chem-viewer-grid' : 'chem-diversity-search';
    const el = document.getElementsByClassName(className);
    const errEl = document.getElementsByClassName('chem-malformed-molecule-error');
    return (el.length && el[0].children.length === 12) || !!errEl.length;
  }, 'Viewer hasn\'t been rendered', 3000);
}

async function awaitScaffoldTree(viewer: DG.Viewer): Promise<void> {
  await awaitCheck(() => document.getElementsByClassName('mol-host').length > 0,
    'scaffold tree has not been generated', 150000);
}

const awaitViewers: {[key: string]: (viewer: DG.Viewer) => Promise<void>} = {
  'Chem Similarity Search': awaitChemSearchViewer,
  'Chem Diversity Search': awaitChemSearchViewer,
  'Scaffold Tree': awaitScaffoldTree
}


category('viewers', () => {
  const df = grok.data.demo.molecules(15);
  const viewers = DG.Func.find({ package: 'Chem', tags: ['viewer'] }).reduce<string[]>((result, f) => {
    const name = f.friendlyName;
    if (name !== 'Scaffold Tree')
      result.push(name);
    return result;
  }, []);

  for (const v of viewers) {
    test(v, async () => {//@ts-ignore
      await testViewer(v, df.clone(), {
        detectSemanticTypes: true,
        readOnly: false,
        arbitraryDfTest: false,
        awaitViewer: awaitViewers[v]
      });
    }, {timeout: 90000});
  }
});