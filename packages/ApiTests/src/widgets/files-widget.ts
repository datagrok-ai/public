// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';
// import $ from 'cash-dom';
// import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';


// category('Widgets', () => {
//   let testConnection: DG.DataConnection;
//   let packageDataConnection: DG.DataConnection;
//   const timeout = 2500;
//   const labelSelector = '.d4-tree-view-tri.d4-tree-view-tri-expanded + i + .d4-tree-view-group-label';

//   before(async () => {
//     testConnection = await grok.dapi.connections.filter('shortName = "Home"').first();
//     packageDataConnection = await grok.dapi.connections.filter('shortName = "AppData"').first();
//   });

//   test('Files', async () => {
//     const fw = DG.FilesWidget.create();
//     expect(fw instanceof DG.FilesWidget, true);
//     expect(fw.root instanceof HTMLElement, true);
//     expect(fw.root.classList.contains('d4-tree-view-root'), true);
//     expect(ui.fileBrowser() instanceof DG.FilesWidget, true)

//     if (testConnection) {
//       const testFW = ui.fileBrowser({path: testConnection.nqName});
//       setTimeout(() => {
//         const label = $(testFW.root).find(labelSelector)[0];
//         expect(label != null, true);
//         expect(label!.textContent, testConnection.friendlyName);
//       }, timeout);
//     }

//     if (packageDataConnection) {
//       const packageName = 'ApiTests';
//       const packageDir = 'datasets';
//       const testFW = ui.fileBrowser({path: `${packageDataConnection.nqName}/${packageName}/${packageDir}`});
//       setTimeout(() => {
//         const labels = $(testFW.root).find(labelSelector);
//         expect(labels[0] != null, true);
//         expect(labels[0]!.textContent, packageDataConnection.friendlyName);
//         expect(labels[1] != null, true);
//         expect(labels[1]!.textContent, packageName);
//         expect(labels[2] != null, true);
//         expect(labels[2]!.textContent, packageDir);
//       }, timeout);
//     }
//   });
// });
