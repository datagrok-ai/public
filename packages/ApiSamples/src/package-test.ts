import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { runTests, tests, TestContext, category, test as _test, delay, initAutoTests as initCoreTests, expect, awaitCheck, before } from '@datagrok-libraries/utils/src/test';
export const _package = new DG.Package();
export { tests };

import './tests/test';

const skip = [
  // Skipped
  'function-events', 'demo', 'ui-events', 'last-error', 'chem-benchmark',
  'menu-customization', '10k-columns-updates', '100-million-rows',
  'files' /* do not test manually */,

  // To fix
  'custom-viewer-properties',
  'open-table-by-id',
  'charts-in-cells',
  'property-grid',
  'tree-view-adv',
  'attached-properties',
  'all-input-types',
  'add-single-filter',
  'custom-filters',
  'filter-group',
  'dynamic-loading',
  //performance 
  'read-strings',
  'dataframe-access',
  '1m-aggregation',
  '100-million-rows',
  '1-million-columns',
  'network-diagram'
];

const scriptViewer = [
  'parameterValidation',
  'parameter expressions',
  'docking', 
  'input-api',
  'helm-input-ui',
  'output-layouts'
];

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  testContext = new TestContext(false, false);
  const data = await runTests({ category, test, testContext });
  return DG.DataFrame.fromObjects(data)!;
}
interface ScriptObject {
  [key: string]: () => Promise<void>;
}

let beforeArr : ScriptObject= {
  ['Scripts:ui:inputs'] : async () => {await grok.functions.call("Helm:getHelmHelper()");}
}

//tags: init
export async function initTests() {

  const scripts = await grok.dapi.scripts.filter('package.shortName = "ApiSamples"').list();
  for (const script of scripts) {
    let catName = ('Scripts:' + script.options.path as string).replaceAll('/', ':');
    category(catName, () => {
      
      if(beforeArr[catName]){
        let currentBefore = beforeArr[catName];
        before(async ()=>{
          await currentBefore();
        })
        delete beforeArr[catName];
      }
      
      _test(script.friendlyName, async () => {

        const annotation = `//name: ${script.friendlyName} 
        //language: javascript
        `;
        debugger
        if (scriptViewer.includes(script.friendlyName))
          await runScriptViewer(script);
        else
          await evaluateScript(script);
        grok.shell.closeAll();

        async function runScriptViewer(script: DG.Script) { 
          const scriptResult = new Promise<boolean>(async (resolve) => { 
            script.script = `${annotation}\n${script.script}`;
            const scriptView = DG.ScriptView.create(script);
            grok.shell.addView(scriptView);
            let timeout: any;
            const subscription = grok.functions.onAfterRunAction.subscribe((funcCall) => {
              if ((funcCall.func as any).script === script.script.replaceAll('\r', '')) {
                if (timeout)
                  clearTimeout(timeout);
                subscription.unsubscribe();
                scriptView.close();
                resolve(true);
              }
            });
            await delay(1000);
            (document.getElementsByClassName("fa-play")[0] as any).click();
            await delay(1000);
            timeout = setTimeout(() => {
              subscription.unsubscribe();
              scriptView.close();
              resolve(false);
            }, 10000);
          })

          if (!(await scriptResult))
            throw new Error(`Script ${'Scripts:' + script.options.path as string}.${script.friendlyName}`);
        }

        async function evaluateScript(script: DG.Script) {
          let timeout: any;
          const subscription = grok.functions.onAfterRunAction.subscribe((funcCall) => {
            if ((funcCall.func as any).script === script.script.replaceAll('\r', '')) {
              if (timeout)
                clearTimeout(timeout);
              subscription.unsubscribe();  
            }
          });
          await script.apply();
          timeout = setTimeout(() => {
            subscription.unsubscribe();
            throw new Error('Script didnt pass');
          }, 10000);
        }
      }, skip.includes(script.friendlyName) ? { skipReason: 'skip' } : { timeout: 60000 });
    });
  }
}

//name: initAutoTests
export async function initAutoTests() {
  await initCoreTests(_package, _package.getModule('package-test.js'));
}
