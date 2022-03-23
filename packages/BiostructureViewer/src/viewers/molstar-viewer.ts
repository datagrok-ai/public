import * as grok from 'datagrok-api/grok';

//@ts-ignore
import {Viewer} from 'molstar/build/viewer/molstar';

// TODO: find out which extensions are needed.
/*const Extensions = {
    'cellpack': PluginSpec.Behavior(CellPack),
    'dnatco-confal-pyramids': PluginSpec.Behavior(DnatcoConfalPyramids),
    'pdbe-structure-quality-report': PluginSpec.Behavior(PDBeStructureQualityReport),
    'rcsb-assembly-symmetry': PluginSpec.Behavior(RCSBAssemblySymmetry),
    'rcsb-validation-report': PluginSpec.Behavior(RCSBValidationReport),
    'anvil-membrane-orientation': PluginSpec.Behavior(ANVILMembraneOrientation),
    'g3d': PluginSpec.Behavior(G3DFormat),
    'model-export': PluginSpec.Behavior(ModelExport),
    'mp4-export': PluginSpec.Behavior(Mp4Export),
    'geo-export': PluginSpec.Behavior(GeometryExport),
    'ma-quality-assessment': PluginSpec.Behavior(MAQualityAssessment),
};*/

const DefaultViewerOptions = {
  extensions: [],
  layoutIsExpanded: false,
  layoutShowControls: false,
  layoutShowRemoteState: false,
  layoutShowSequence: false,
  layoutShowLog: false,
  layoutShowLeftPanel: false,
  collapseLeftPanel: true,
  collapseRightPanel: true,

  viewportShowExpand: false,
  viewportShowControls: false,
};

/**
 * Creates an instance of Mol* viewer.
 *
 * @param {string} pdbID ID of an entity in PDB.
 */
export async function createViewer(pdbID: string) {
  const view = grok.shell.newView('Mol*');
  const viewerContainer = view.root;
  const v = new Viewer(viewerContainer, DefaultViewerOptions);

  await v.loadPdb(pdbID);
  //v.handleResize();
}

