// //name: PDB Viewer
// //description: 3D structure data for large biological molecules (proteins, DNA, and RNA)
// //top-menu: Bio | PDB ...
// //output: viewer result
// export async function pdbViewer(): Promise<void> {
//   const view: DG.TableView = grok.shell.tv;
//   const pdbTag = view.dataFrame.getTag(pdbTAGS.PDB);
//   if (pdbTag) {
//     const viewer = (await view.dataFrame.plot.fromType('Ngl', {})) as DG.JsViewer;
//     view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL viewer', 0.4);
//   } else {
//     const bsView: DG.TableView = grok.shell.tv;
//
//     const ligandSelection = {};
//     let ligands = ['R', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'S', 'T', 'U', 'V', 'W'];
//     for (let i = 0; i < ligands.length; i++)
//       ligandSelection[ligands[i]] = [false, 400 + i];
//
//     return new Promise<void>((resolve, reject) => {
//       const fileBrowser = ui.fileBrowser({path: `System:AppData/${_package.name}/samples`});
//       const dlg: DG.Dialog = ui.dialog({title: 'Open PDB file'})
//         .add(fileBrowser.root)
//         .addButton('OK', () => {
//           setTimeout(async () => {
//             const filePath: string = fileBrowser.props.file;
//             const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('PDB Viewer');
//             try {
//               const pdbStr: string = await grok.dapi.files.readAsText(filePath);
//
//               const viewer: DG.JsViewer = (await view.dataFrame.plot.fromType('Ngl',
//                 {[nglPROPS.pdb]: pdbStr})) as DG.JsViewer;
//               view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL viewer', 0.4);
//
//               resolve();
//             } catch (err: any) {
//               const errMsg: string = errorToConsole(err);
//               console.error(errMsg);
//               reject(err.toString());
//             } finally {
//               pi.close();
//               dlg.close();
//             }
//           }, 0 /* next event cycle */);
//         })
//         .show();
//     });
//   }
// }
