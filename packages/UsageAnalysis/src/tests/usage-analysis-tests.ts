import {category, expect, delay, test, before} from "@datagrok-libraries/utils/src/test";
import * as grok from "datagrok-api/grok";

category('usageAnalysis', () => {
  const allViewersToView: {[index: string] : string[]} = {
    'Overview': ['Groups', 'Users', 'Unique Users', 'Packages'],
    'Events': ['Events', 'Packages', 'Package Functions', 'Functions', 'Sources'],
    'Errors': ['Errors', 'Errors', 'Disabled Errors', 'Error Sources'],
    'Function Errors': ['Function Errors', 'Function Errors', 'Function Disabled Errors', 'Packages By Errors'],
    'Users': ['Users', 'Usage', 'Unique Users'],
    'Data': ['Queries', 'Queries', 'Connections', 'Data Sources']
  }

  function changeView(viewName: string) {
    let view = grok.shell.view(viewName);
    if (view === null)
      throw `Can't find view ${viewName}`;

    grok.shell.v = view;
  }

  before(async () => {
    await grok.functions.call("UsageAnalysis:usageAnalysisApp");
  });

  test('openApp', async () => {
    await delay(5000);
    expect(grok.shell.v.name === 'Overview', true);
  });

  test('viewsTest', async () => {
    expect(Object.keys(allViewersToView).every(viewName => grok.shell.view(viewName) !== undefined), true);
  });

  test('allViewersTest', async () => {
    for (let viewName of Object.keys(allViewersToView)) {
      changeView(viewName);

      let foundedViewersOfView = 0;

      let h1Elements = document.getElementsByTagName('h1');
      for (let heI = 0; heI < h1Elements.length; heI++) {
        for (let vovI = 0; vovI < allViewersToView[viewName].length; vovI++) {
          if (h1Elements[heI].innerText === allViewersToView[viewName][vovI]) {
            foundedViewersOfView++;
            break;
          }
        }
      }
      expect(foundedViewersOfView, allViewersToView[viewName].length);
    }
  });

});