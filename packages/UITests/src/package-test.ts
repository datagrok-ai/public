import './ui/inputs';
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
import './gui/dialogs';
import './gui/files';
import './gui/grid';
import './gui/data-science';
import './gui/project-upload';
// import './gui/viewers/bar-chart';
// import './gui/viewers/box-plot';
// import './gui/viewers/density-plot';
// import './gui/viewers/form';
// import './gui/viewers/histogram';
// import './gui/viewers/leaflet';
// import './gui/viewers/line-chart';
// import './gui/viewers/matrix-plot';
// import './gui/viewers/network-diagram';
// import './gui/viewers/pc-plot';
// import './gui/viewers/pie-chart';
// import './gui/viewers/scatter-plot';
// import './gui/viewers/scatter-plot-3d';
// import './gui/viewers/word-cloud';

import * as DG from 'datagrok-api/dg';
import {runTests, tests} from '@datagrok-libraries/utils/src/test';

export const _package = new DG.Package();
export {tests};

//name: test
//output: dataframe result
export async function test(): Promise<DG.DataFrame> {
  const data = await runTests();
  return DG.DataFrame.fromObjects(data)!;
}
