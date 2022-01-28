import {category, expect, delay, test} from "@datagrok-libraries/utils/src/test";
import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";

category('usageAnalysis', () => {
  const allViewersToView: {[index: string] : string[]} = {
    'Overview': ['Total Users', 'Unique Users List', 'Unique Users', 'Events', 'Errors'],
    'Events': ['Events', 'Packages', 'Package Functions', 'Functions', 'Sources'],
    'Errors': ['Errors', 'Errors', 'Error Sources'],
    'Users': ['Users', 'Usage', 'Unique Users'],
    'Data': ['Queries', 'Queries', 'Connections', 'Data Sources']
  }

  const canvasViewersToView: {[index: string] : {[index: string] : number}} = {
    'Overview': {'d4-line-chart': 3},
    'Events': {'d4-line-chart': 1, 'd4-bar-chart': 4},
    'Errors': {'d4-line-chart': 1, 'd4-bar-chart': 3},
    'Users': {'d4-line-chart': 1, 'd4-scatter-plot': 1, 'd4-bar-chart': 1},
    'Data': {'d4-line-chart': 1, 'd4-bar-chart': 3}
  }

  function changeView(viewName: string) {
    let view = grok.shell.view(viewName);
    if (view === null)
      throw `Can't find view ${viewName}`;

    grok.shell.v = view;
  }

  test('openApp', async () => {
    await grok.functions.call("UsageAnalysis:usageAnalysisApp");
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

  test('canvasViewersTest', async () => {
    for (let viewName of Object.keys(canvasViewersToView)) {
      changeView(viewName);
      await delay(2000);

      for (let viewerName of Object.keys(canvasViewersToView[viewName])) {
        if (document.getElementsByClassName(viewerName).length !== canvasViewersToView[viewName][viewerName])
          throw `Not enough ${viewerName} in ${viewName}`;
      }
    }
  });

});