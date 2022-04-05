/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './dataframe/data_frame';
import './dataframe/calculated-columns';
import './dataframe/events';
import './shell/shell';
import './shell/windows';
import './viewer/viewer';
import './views/docking';
import './views/events';
import './views/layouts';
import './dapi/files';
import './dapi/fetch';
import './dapi/admin';
import './dapi/groups';
import './ui/inputs';
import './ui/forms';
import './dapi/dapi';
import './dapi/entities';
import './dapi/layouts';
import './dapi/projects';
import './dapi/tables';
import './dapi/user-data-storage';
import './dapi/users';
import './shell/ml';
import './ui/divs';
import './ui/buttons';
import './widgets/legend';
import './ui/icons';
import './ui/tables';
import './ui/rangeSlider';
import './ui/accordion';
import './ui/tabControl';
import './ui/list';
import './ui/image';
import './ui/viewers-adding';
import './grid/grid';

import {runTests} from "@datagrok-libraries/utils/src/test";
export let _package = new DG.Package();


//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//output: dataframe result
//top-menu: Tools | Dev | JS API Tests
export async function test(category: string, test: string): Promise<DG.DataFrame> {
  let data = await runTests({category, test});
  return DG.DataFrame.fromObjects(data)!;
}

//name: testPackages
//output: dataframe result
//top-menu: Tools | Dev | Test Packages
export async function testPackages(): Promise<DG.DataFrame> {
  let funcs = DG.Func.find({name:'test'});
  let dfs:DG.DataFrame[] = [];
  for (const f of funcs) {
    if (f.package?.name != null) {
      grok.shell.closeAll();
      grok.shell.info(`Testing ${f.package.name}`);
      let df = await f.apply();
      if (df == null) {
        grok.shell.error(`Failed to fetch test results from ${f.package.name}`);
        continue;
      }
      let packageColumn = DG.Column.string('package', df.rowCount);
      packageColumn.init((n) => f.package.name);
      df.columns.insert(packageColumn, 0);
      dfs.push(df);
      grok.shell.closeAll();
    }
  }

  let result: DG.DataFrame | null = null;
  for (const df of dfs) {
    if (result == null)
      result = df;
    else result.append(df, true);
  }

  return result!;
}



//name: testManager
//top-menu: Tools | Dev | Test Manager
export async function testManager() {
  let testFunctions = DG.Func.find({ name: 'Test' });
  const packagesTestsList = {};
  for (let f of testFunctions) {
    //@ts-ignore
    await f.package.load({ file: f.options.file });
    //@ts-ignore
    const allPackageTests = f.package.getModule(f.options.file).tests;
    if (allPackageTests) {
      //@ts-ignore
      packagesTestsList[ f.package.friendlyName ] = allPackageTests;
    }
  }
  createTestMangerUI(packagesTestsList);
}

function createTestMangerUI(packagesTests: any) {
  Object.keys(packagesTests).forEach(pack => {
    Object.keys(packagesTests[ pack ]).forEach(cat => {
      //@ts-ignore
      packagesTests[ pack ][ cat ].tests.forEach(t => t.active = true);
    })
  });
  const v = grok.shell.newView('Test manager');
  const acc = ui.accordion();
  Object.keys(packagesTests).forEach(pack => {
    acc.addCountPane(pack, () => {
      const catAcc = ui.accordion();
      Object.keys(packagesTests[ pack ]).forEach(cat => {
        catAcc.addCountPane(cat, () => {
          let testsDf = DG.DataFrame.create(packagesTests[ pack ][ cat ].tests.length);
          //@ts-ignore
          testsDf.columns.addNewString('Test').init((i) => packagesTests[ pack ][ cat ].tests[ i ].name);
          //@ts-ignore
          testsDf.columns.addNewBool('Active').init((i) => packagesTests[ pack ][ cat ].tests[ i ].active);
          testsDf.onCurrentCellChanged.subscribe(() => {
            setTimeout(() => {
              if (testsDf.currentCol.name === 'Active') {
                packagesTests[ pack ][ cat ].tests[ testsDf.currentRowIdx ].active = testsDf.currentCol.get(testsDf.currentRowIdx);
                //@ts-ignore
                catAcc.getPane(`${cat}`).root.children[ 0 ].children[ 0 ].innerHTML = `${packagesTests[ pack ][ cat ].tests.filter(it => it.active).length}`;
              }
            }, 100);
          });
          return testsDf.plot.grid().root;
        },
          //@ts-ignore
          () => packagesTests[ pack ][ cat ].tests.filter(it => it.active).length);
        let panel = catAcc.getPane(`${cat}`);
        //@ts-ignore
        $(panel.root).css('display', 'flex');
        //@ts-ignore
        $(panel.root).css('opacity', '1');
      });
      return catAcc.root;
    }, () => Object.keys(packagesTests[ pack ]).length);
  });
  const runTestsButton = ui.bigButton('Run', () => {
    Object.keys(packagesTests).forEach(pack => {
      Object.keys(packagesTests[ pack ]).forEach(cat => {
        //@ts-ignore
        packagesTests[ pack ][ cat ].tests.forEach(t => {
          if (t.active) {
            t.test();
          }
        });
      })
    });
  });
  v.append(runTestsButton);
  v.append(acc.root);
}

