import * as DG from 'datagrok-api/dg';
import {runTests, tests, TestContext, initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';

import './ui/inputs';
import './ui/input-for-property';
import './ui/inputform';
import './ui/forms';
import './ui/divs';
import './ui/buttons';
import './ui/icons';
import './ui/tables';
import './ui/range-slider';
import './ui/accordion';
import './ui/tab-control';
import './ui/list';
import './ui/image';
import './ui/users';
import './ui/groups';
import './ui/tags';
import './ui/sharing';
import './ui/get-all-top100';
import './gui/dialogs';
import './gui/files';
import './gui/grid';
import './gui/project-upload';
import './views/docking';
import './views/docking-nested';
import './views/events';
import './views/layouts';
import './views/files-view';
import './viewers/viewers';
import './viewers/filters';
import './gui/viewers/scatter-plot';
import './shell/windows';

// import './gui/apps';

// import './gui/viewers/bar-chart';
// import './gui/viewers/box-plot';
// import './gui/viewers/density-plot';
import './gui/viewers/form';
import './gui/viewers/histogram';
// import './gui/viewers/leaflet';
import './gui/viewers/line-chart';
// import './gui/viewers/matrix-plot';
// import './gui/viewers/network-diagram';
// import './gui/viewers/pc-plot';
// import './gui/viewers/pie-chart';
// import './gui/viewers/scatter-plot-3d';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  testContext = new TestContext(false, false);
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
