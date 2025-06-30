import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:Init', {});
  }

  export async function molecule3dCellRenderer(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:Molecule3dCellRenderer', {});
  }

  export async function pdbIdCellRenderer(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PdbIdCellRenderer', {});
  }

  export async function viewPdbById(pdbId: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ViewPdbById', { pdbId });
  }

  export async function viewPdbByData(pdbData: string, name: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ViewPdbByData', { pdbData, name });
  }

  export async function getNglGlService(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:GetNglGlService', {});
  }

  export async function viewBiostructure(content: string, format: string, name: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ViewBiostructure', { content, format, name });
  }

  //Opens PDB file
  export async function importPdb(fileContent: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ImportPdb', { fileContent });
  }

  //Opens XYZ file
  export async function importXYZ(fileContent: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ImportXYZ', { fileContent });
  }

  //Opens biostructure files supported with NGL
  export async function importWithNgl(fileContent: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ImportWithNgl', { fileContent });
  }

  //Opens .pdbqt file with docking result ligand poses
  export async function importPdbqt(fileContent: string, test: boolean): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ImportPdbqt', { fileContent, test });
  }

  export async function previewNglStructure(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PreviewNglStructure', { file });
  }

  export async function previewNglSurface(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PreviewNglSurface', { file });
  }

  export async function previewNglDensity(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PreviewNglDensity', { file });
  }

  export async function previewBiostructureStructure(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PreviewBiostructureStructure', { file });
  }

  export async function previewBiostructureTopology(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PreviewBiostructureTopology', { file });
  }

  export async function previewBiostructureDensity(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PreviewBiostructureDensity', { file });
  }

  export async function openPdbResidues(fi: DG.FileInfo): Promise<any> {
    return await grok.functions.call('BiostructureViewer:OpenPdbResidues', { fi });
  }

  export async function pdbIdNglPanelWidget(pdbId: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:PdbIdNglPanelWidget', { pdbId });
  }

  //Example app for NGL drawing in grid cells
  export async function nglForGridTestApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:NglForGridTestApp', {});
  }

  //Test app for NglViewer
  export async function nglViewerApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:NglViewerApp', {});
  }

  //Test app for BiostructureViewer (molstar)
  export async function biostructureViewerApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:BiostructureViewerApp', {});
  }

  //Test app for BiotrackViewer (saguaro)
  export async function biotrackViewerApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:BiotrackViewerApp', {});
  }

  //Test app for twin BiostructureViewer (molstar) and BiotrackViewer (saguaro)
  export async function biostructureAndTrackViewerApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:BiostructureAndTrackViewerApp', {});
  }

  export async function ligandsWithNglApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:LigandsWithNglApp', {});
  }

  export async function ligandsWithBiostructureApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:LigandsWithBiostructureApp', {});
  }

  export async function biostructureDataProviderApp(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:BiostructureDataProviderApp', {});
  }

  //3D structure viewer for large biological molecules (proteins, DNA, and RNA)
  export async function nglViewer(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:NglViewer', {});
  }

  //3D structure molstar RCSB viewer for large biological molecules (proteins, DNA, and RNA)
  export async function molstarViewer(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:MolstarViewer', {});
  }

  //structure polymer annotation tracks
  export async function saguaroViewer(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:SaguaroViewer', {});
  }

  export async function getPdbHelper(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:GetPdbHelper', {});
  }

  export async function dockingDemo(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:DockingDemo', {});
  }

  export async function inGridDemo(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:InGridDemo', {});
  }

  export async function copyRawBiostructureValue(gridCell: any): Promise<any> {
    return await grok.functions.call('BiostructureViewer:CopyRawBiostructureValue', { gridCell });
  }

  export async function downloadRawBiostructureValue(gridCell: any): Promise<any> {
    return await grok.functions.call('BiostructureViewer:DownloadRawBiostructureValue', { gridCell });
  }

  export async function showBiostructureViewerMenuItem(gridCell: any): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ShowBiostructureViewerMenuItem', { gridCell });
  }

  export async function showNglViewerMenuItem(gridCell: any): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ShowNglViewerMenuItem', { gridCell });
  }

  export async function openTableResiduesMenuItem(fi: any): Promise<any> {
    return await grok.functions.call('BiostructureViewer:OpenTableResiduesMenuItem', { fi });
  }

  //Display ligand poses along the structure
  export async function demoBioDockingConformations(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:DemoBioDockingConformations', {});
  }

  //View structures PDB in grids
  export async function demoBioProteins(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:DemoBioProteins', {});
  }

  export async function demoFix1(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:DemoFix1', {});
  }

  export async function demoFix2(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:DemoFix2', {});
  }

  export async function demoFix3(): Promise<any> {
    return await grok.functions.call('BiostructureViewer:DemoFix3', {});
  }

  export async function readAsText(file: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ReadAsText', { file });
  }

  export async function readAsTextDapi(file: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:ReadAsTextDapi', { file });
  }

  //Get biostructure by id as PDB
  export async function getBiostructureRcsbPdb(id: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:GetBiostructureRcsbPdb', { id });
  }

  //Get biostructure by id as mmCIF
  export async function getBiostructureRcsbMmcif(id: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:GetBiostructureRcsbMmcif', { id });
  }

  //Get biostructure by id as BinaryCIF
  export async function getBiostructureRcsbBcif(id: string): Promise<any> {
    return await grok.functions.call('BiostructureViewer:GetBiostructureRcsbBcif', { id });
  }

  //Packs BiostructureData value into JSON string
  export async function biostructureDataToJson(binary: boolean, data: any, ext: string, options: any): Promise<any> {
    return await grok.functions.call('BiostructureViewer:BiostructureDataToJson', { binary, data, ext, options });
  }

  export async function structure3D(molecule: any): Promise<any> {
    return await grok.functions.call('BiostructureViewer:Structure3D', { molecule });
  }
}
