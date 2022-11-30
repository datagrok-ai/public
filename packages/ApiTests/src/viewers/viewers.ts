// import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

// import {after, before, category, expect, test, delay} from '@datagrok-libraries/utils/src/test';
// import {TestViewerForProperties} from './test-viewer-for-properties';
// import {waitForElement} from '../gui/gui-utils';


// category('Viewers', () => {
//   let df: DG.DataFrame;
//   let tv: DG.TableView;
//   let coreViewerTypes: string[];

//   let viewerList: DG.JsViewer[];

//   before(async () => {
//     coreViewerTypes = Object.values(DG.VIEWER).filter((v) =>
//       v != DG.VIEWER.TIMELINES &&
//       v != DG.VIEWER.GLOBE &&
//       v != DG.VIEWER.SCATTER_PLOT_3D &&
//       v != DG.VIEWER.GOOGLE_MAP &&
//       v != DG.VIEWER.SHAPE_MAP &&
//     v != DG.VIEWER.SURFACE_PLOT);
//     df = grok.data.demo.demog();
//     tv = grok.shell.addTableView(df);

//     viewerList = [];
//   });

//   test('addViewer(ViewerType)', async () => {
//     try {
//       for (const viewerType of coreViewerTypes) {
//         console.log(`Adding ${viewerType}`);
//         const viewer = await addViewerAndWait(tv, viewerType);//tv.addViewer(viewerType);
//         if (!(viewer instanceof DG.Viewer))
//           throw `TableView.addViewer('${viewerType}') should add a Viewer instance`;
//         await delay(200);

//         viewer.removeFromView();
//         console.log(`Removed ${viewerType}`);
//       }
//     } finally {
//       closeViewers(tv);
//     }
//   });

//   test('addViewer(Viewer)', async () => {
//     try {
//       for (const viewerType of coreViewerTypes) {
//         console.log(`Adding ${viewerType}`);
//         const viewer = await addViewerAndWait(tv, DG.Viewer.fromType(viewerType, df));
//         await delay(200);
//         viewer.removeFromView();
//         console.log(`Removed ${viewerType}`);
//       }
//     } finally {
//       closeViewers(tv);
//     }
//   });

//   test('close', async () => {
//     tv.scatterPlot();
//     tv.barChart();
//     expect(Array.from(tv.viewers).length, 3);
//     closeViewers(tv);
//     expect(Array.from(tv.viewers).length, 1);
//   });

//   test('getViewerTypes', async () => {
//     const registeredViewers = DG.Viewer.getViewerTypes();
//     expect(coreViewerTypes.every((t) => registeredViewers.includes(t)), true);
//   });

//   test('fromType', async () => {
//     for (const viewerType of coreViewerTypes) {
//       const viewer = DG.Viewer.fromType(viewerType, df);
//       if (!(viewer instanceof DG.Viewer))
//         throw `Viewer.fromType('${viewerType}', df) should add a Viewer instance`;
//       expect(viewer.table.id, df.id);
//     }
//   });

//   test('ScatterPlotViewer.zoom', async () => {
//     const sp = DG.Viewer.scatterPlot(df);
//     tv.addViewer(sp);
//     try {
//       sp.zoom(10, 10, 100, 100);
//     } finally {
//       sp.close();
//     }
//   });

//   test('ScatterPlotViewer.onZoomed', async () => {
//     const sp = DG.Viewer.scatterPlot(df);
//     tv.addViewer(sp);
//     try {
//       let rectangle: any;
//       sp.onZoomed.subscribe((r) => rectangle = r);
//       sp.zoom(10, 10, 100, 100);
//       expect(rectangle instanceof DG.Rect, true);
//     } finally {
//       sp.close();
//     }
//   });

//   test('onTableAttached', async () => {
//     const df1 = DG.DataFrame.fromCsv('id1');
//     const df2 = DG.DataFrame.fromCsv('id2');

//     const viewer: TestViewerForProperties = (await df1.plot.fromType(
//       'TestViewerForProperties', {})) as TestViewerForProperties;
//     expect(viewer.dataFrame.columns.byIndex(0).name, df1.columns.byIndex(0).name);
//     expect(viewer.onTableAttachedCounter, 1);

//     viewer.dataFrame = df2;
//     expect(viewer.dataFrame.columns.byIndex(0).name, df2.columns.byIndex(0).name);
//     expect(viewer.onTableAttachedCounter, 2);

//     // TODO: Check onTableAttached has been called
//   }, {skipReason: 'GROK-11484'});

//   test('setPropertyStringWithNumber', async () => {
//     // const v: TestViewerForProperties = tv.addViewer('TestViewerForProperties', {}) as TestViewerForProperties;
//     const dfTemp = DG.DataFrame.fromCsv('id');
//     const viewer: TestViewerForProperties = (
//       await dfTemp.plot.fromType('TestViewerForProperties')) as TestViewerForProperties;
//     viewerList.push(viewer);

//     let exCaught: boolean = false;
//     try {
//       viewer.setOptions({'testPropertyString': 1});
//     } catch {
//       // There should be an exception caught while
//       // assigning number value to string option of JsViewer
//       exCaught = true;
//     }
//     const propValueFromProps = viewer.props.get('testPropertyString');
//     const propValueFromObject = viewer.testPropertyString;

//     // Silent type cast or exception expected
//     if ((typeof propValueFromObject !== typeof 'str' || typeof propValueFromProps !== typeof 'str') && !exCaught)
//       throw new Error('JsViewer string property assigned with number value ' +
//         `become value of type '${typeof propValueFromObject}' without an exception or type conversion.`);
//   }, {skipReason: 'GROK-11485'});

//   test('setPropertyIntWithString', async () => {
//     // const v: TestViewerForProperties = tv.addViewer('TestViewerForProperties', {}) as TestViewerForProperties;
//     const dfTemp = DG.DataFrame.fromCsv('id');
//     const viewer: TestViewerForProperties = (
//       await dfTemp.plot.fromType('TestViewerForProperties')) as TestViewerForProperties;
//     viewerList.push(viewer);

//     let exCaught: boolean = false;
//     try {
//       viewer.setOptions({'testPropertyInt': '1'}); // silent parse to int available
//     } catch {
//       // There should be an exception caught while
//       // assigning string value to int option of JsViewer
//       exCaught = true;
//     }
//     const propValueFromProps = viewer.props.get('testPropertyInt');
//     const propValueFromObject = viewer.testPropertyInt;

//     // Silent type cast or exception expected
//     if ((typeof propValueFromObject !== typeof 1 || typeof propValueFromProps !== typeof 1) && !exCaught)
//       throw new Error('JsViewer int property assigned with string value ' +
//         `become value of type '${typeof propValueFromObject}' without an exception or type conversion.`);
//   }, {skipReason: 'GROK-11485'});


//   test('testViewersLayout', async () => {
//     let viewers = DG.Viewer.getViewerTypes().filter(vt => !vt.toLowerCase().includes('widget'));
//     const skipViewers = ['3d scatter plot', 'Google map', 'Network diagram', 'Sankey'];
//     let layout, res, errorViewers = [], i = 0;
    
//     for (let v of viewers) {
//       if (skipViewers.includes(v)) continue;
//       res = [];
//       try { 
//         tv.addViewer(v);
//         await delay(100);
//         res.push(Array.from(tv.viewers).length);
//         layout = tv.saveLayout();
//         tv.resetLayout();
//         res.push(Array.from(tv.viewers).length);
//         tv.loadLayout(layout);
//         await delay(1000);
//         res.push(Array.from(tv.viewers).length);
//         await testBoolProps();
//       } catch (e: any) {
//         errorViewers.push([v, e.message]);
//       } finally {
//         i++;
//         console.log(v, i);
//         if (!(res[0] === 2 && res[1] === 1 && res[2] === 2)) errorViewers.push([v, res]);
//         tv.resetLayout();
//       }
//     }

//     grok.shell.info(`Tested ${i} viewers of ${viewers.length}, skipped ${skipViewers.length}, error in ${errorViewers.length}`);
//     if (errorViewers.length !== 0) throw `Error viewers: ${errorViewers}`;
//   }, {skipReason: 'too long execution'});


//   after(async () => {
//     tv.close();
//     grok.shell.closeTable(df);

//     for (const viewer of viewerList) {
//       try {
//         // viewer.removeFromView();
//         viewer.close();
//       } catch (err: any) {
//         console.warn(`Closing viewer error: ${err.toString()}`);
//       } // ignore everything on closing viewers
//     }
//     viewerList = [];
//   });
// });

// function closeViewers(view: DG.TableView) {
//   Array.from(view.viewers).slice(1).forEach((v) => v.close());
// }

// function addViewerAndWait(tv: DG.TableView, viewerType: string | DG.Viewer): Promise<DG.Viewer> {
//   return new Promise((resolve, reject) => {
//     const sub = grok.events.onViewerAdded.subscribe((data) => {
//       // @ts-ignore
//       if ((data.args.viewer as DG.Viewer).type == viewerType || data.args.viewer == viewerType) {
//         sub.unsubscribe();
//         // @ts-ignore
//         resolve(data.args.viewer);
//       }
//     });
//     tv.addViewer(viewerType);
//     setTimeout(() => {
//       sub.unsubscribe();
//       // eslint-disable-next-line prefer-promise-reject-errors
//       reject('timeout');
//     }, 100);
//   });
// }

// async function testBoolProps() {
//   const params = document.querySelectorAll('.grok-font-icon-settings')[2] as HTMLElement;
//   params.click();
//   await delay(1000);
//   const checkboxes = Array.from(document.querySelectorAll('.property-grid-item-editor-checkbox'));
//   for (let cb of checkboxes) {
//     (cb as HTMLElement).click();
//     await delay(50);
//   }
// }
