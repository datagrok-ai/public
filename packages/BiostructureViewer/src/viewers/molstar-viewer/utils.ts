import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Viewer as RcsbViewer, ViewerProps as RcsbViewerProps} from '@rcsb/rcsb-molstar/build/src/viewer';
import {BuiltInTrajectoryFormat, BuiltInTrajectoryFormats} from 'molstar/lib/mol-plugin-state/formats/trajectory';
import {BuildInStructureFormat} from 'molstar/lib/mol-plugin-state/formats/structure';
import {BuildInShapeFormat} from 'molstar/lib/mol-plugin-state/formats/shape';
import {BuildInVolumeFormat} from 'molstar/lib/mol-plugin-state/formats/volume';
import {Unsubscribable} from 'rxjs';
import {BiostructureProps, BiostructurePropsDefault} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {Color as msColor} from 'molstar/lib/mol-util/color';

import {defaults, molecule3dFileExtensions} from './consts';
import {parseAndVisualsData} from './molstar-viewer-open';

export async function initViewer(viewName: string = 'Mol*'): Promise<RcsbViewer> {
  const view = grok.shell.newView(viewName);
  const viewerContainer = view.root;
  const viewer = new RcsbViewer(viewerContainer, castProps(defaults));
  return viewer;
}

/**
 * Creates an instance of Mol* viewer.
 *
 * @param {string} pdbId ID of an entity in PDB.
 */
export async function byId(pdbId: string) {
  await initViewer()
    .then((v: RcsbViewer) => {
      v.loadPdbId(pdbId);
    });
  //v.handleResize();
}

/**
 * Creates an instance of Mol* viewer.
 *
 * @param {string} data Data in PDB
 */
export async function byData(data: string, name: string = 'Mol*', format: BuiltInTrajectoryFormat = 'pdb') {
  await initViewer(name)
    .then(async (viewer: RcsbViewer) => {
      // detecting format by data content
      let binary: boolean = false;
      if (data.includes('mmcif_pdbx.dic')) {
        format = 'mmcif';
        binary = false;
      }

      await viewer.loadStructureFromData(data, format, binary);
    });
  //v.handleResize();
}

export async function viewMolstarUI(content: string, name?: string, format?: BuiltInTrajectoryFormat): Promise<void> {
  await byData(content, name, format);
}

/** Creates view with Molstar viewer to preview Biostructure (PDB)
 * returns the view immediately, but the viewer is created asynchronously and promise
 * for that is returned separately which resolves once the viewer is initialized.
 */
export function previewMolstarUI(file: DG.FileInfo): { view: DG.View, loadingPromise: Promise<void> } {
  const builtinFormats = BuiltInTrajectoryFormats.map((obj) => obj[0]) as string[];

  function isSupportedFormat(fi: DG.FileInfo) {
    return fi.extension in molecule3dFileExtensions;
  }

  if (!isSupportedFormat(file)) {
    grok.shell.error(`Unsupported format: ${file.extension}`);
    throw new Error(`Unsupported format: ${file.extension}`);
  }

  const view = DG.View.create();
  const viewer = new RcsbViewer(view.root, castProps(defaults));
  const subs: Unsubscribable[] = [];
  subs.push(ui.onSizeChanged(view.root).subscribe((value: any) => {
    viewer.handleResize();
  }));
  subs.push(grok.events.onViewRemoved.subscribe((evtView) => {
    if (evtView.id === view.id)
      for (const sub of subs) sub.unsubscribe();
  }));

  const loadingPromise = new Promise<void>(async (resolve, reject) => {
    try {
      const { binary} = molecule3dFileExtensions[file.extension];
      const data: string | Uint8Array = binary ? await file.readAsBytes() : await file.readAsString();
      await parseAndVisualsData(viewer.plugin, {ext: file.extension, data: data});
      resolve();
    } catch (err: any) {
      reject(err);
    }
  });

  return {view, loadingPromise};
}

function castProps(src: BiostructureProps): Partial<RcsbViewerProps> {
  const res: Partial<RcsbViewerProps> = {
    showImportControls: src.showImportControls,
    showWelcomeToast: src.showWelcomeToast,
    layoutIsExpanded: src.layoutIsExpanded,
    layoutShowLog: src.layoutShowLog,
    layoutControlsDisplay: src.layoutControlsDisplay,
    layoutShowControls: src.layoutShowControls,
    layoutShowSequence: src.layoutShowSequence,
    backgroundColor: src.backgroundColor as msColor,
    detachedFromSierra: src.detachedFromSierra,
    showExportControls: src.showExportControls,
    showMembraneOrientationPreset: src.showMembraneOrientationPreset,
    showSessionControls: src.showSessionControls,
    showStructureSourceControls: src.showStructureSourceControls,
    viewportShowExpand: src.viewportShowExpand,
    showSuperpositionControls: src.showSuperpositionControls,
    viewportShowSelectionMode: src.viewportShowSelectionMode,
    volumeStreamingServer: src.volumeStreamingServer,
    modelUrlProviders: [],
    extensions: []
  };
  return res;
}
