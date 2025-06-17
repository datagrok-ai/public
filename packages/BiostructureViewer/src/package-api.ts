import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  //export async function init() 
  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:Init', {});
  }

  export async function molecule3dCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:Molecule3dCellRenderer', {});
  }

  export async function pdbIdCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PdbIdCellRenderer', {});
  }

  export async function viewPdbById(pdbId: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ViewPdbById', { pdbId });
  }

  export async function viewPdbByData(pdbData: string, name: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ViewPdbByData', { pdbData, name });
  }

  export async function getNglGlService(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:GetNglGlService', {});
  }

  export async function viewBiostructure(content: string, format: string, name: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ViewBiostructure', { content, format, name });
  }

  //Opens PDB file
  export async function importPdb(fileContent: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ImportPdb', { fileContent });
  }

  //Opens XYZ file
  export async function importXYZ(fileContent: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ImportXYZ', { fileContent });
  }

  //Opens biostructure files supported with NGL
  export async function importWithNgl(fileContent: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ImportWithNgl', { fileContent });
  }

  //Opens .pdbqt file with docking result ligand poses
  export async function importPdbqt(fileContent: string, test: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ImportPdbqt', { fileContent, test });
  }

  export async function previewNglStructure(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PreviewNglStructure', { file });
  }

  export async function previewNglSurface(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PreviewNglSurface', { file });
  }

  export async function previewNglDensity(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PreviewNglDensity', { file });
  }

  export async function previewBiostructureStructure(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PreviewBiostructureStructure', { file });
  }

  export async function previewBiostructureTopology(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PreviewBiostructureTopology', { file });
  }

  export async function previewBiostructureDensity(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PreviewBiostructureDensity', { file });
  }

  export async function openPdbResidues(fi: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:OpenPdbResidues', { fi });
  }

  export async function pdbIdNglPanelWidget(pdbId: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:PdbIdNglPanelWidget', { pdbId });
  }

  //Example app for NGL drawing in grid cells
  export async function nglForGridTestApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:NglForGridTestApp', {});
  }

  //Test app for NglViewer
  export async function nglViewerApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:NglViewerApp', {});
  }

  //Test app for BiostructureViewer (molstar)
  export async function biostructureViewerApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:BiostructureViewerApp', {});
  }

  //Test app for BiotrackViewer (saguaro)
  export async function biotrackViewerApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:BiotrackViewerApp', {});
  }

  //Test app for twin BiostructureViewer (molstar) and BiotrackViewer (saguaro)
  export async function biostructureAndTrackViewerApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:BiostructureAndTrackViewerApp', {});
  }

  export async function ligandsWithNglApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:LigandsWithNglApp', {});
  }

  export async function ligandsWithBiostructureApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:LigandsWithBiostructureApp', {});
  }

  export async function biostructureDataProviderApp(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:BiostructureDataProviderApp', {});
  }

  //3D structure viewer for large biological molecules (proteins, DNA, and RNA)
  export async function nglViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:NglViewer', {});
  }

  //3D structure molstar RCSB viewer for large biological molecules (proteins, DNA, and RNA)
  export async function molstarViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:MolstarViewer', {});
  }

  //structure polymer annotation tracks
  export async function saguaroViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:SaguaroViewer', {});
  }

  export async function getPdbHelper(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:GetPdbHelper', {});
  }

  //export async function dockingDemo() 
  export async function dockingDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:DockingDemo', {});
  }

  //export async function inGridDemo() 
  export async function inGridDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:InGridDemo', {});
  }

  export async function copyRawBiostructureValue(gridCell: any): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:CopyRawBiostructureValue', { gridCell });
  }

  export async function downloadRawBiostructureValue(gridCell: any): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:DownloadRawBiostructureValue', { gridCell });
  }

  export async function showBiostructureViewerMenuItem(gridCell: any): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ShowBiostructureViewerMenuItem', { gridCell });
  }

  export async function showNglViewerMenuItem(gridCell: any): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ShowNglViewerMenuItem', { gridCell });
  }

  export async function openTableResiduesMenuItem(fi: any): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:OpenTableResiduesMenuItem', { fi });
  }

  //Display ligand poses along the structure
  export async function demoBioDockingConformations(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:DemoBioDockingConformations', {});
  }

  //View structures PDB in grids
  export async function demoBioProteins(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:DemoBioProteins', {});
  }

  export async function demoFix1(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:DemoFix1', {});
  }

  export async function demoFix2(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:DemoFix2', {});
  }

  export async function demoFix3(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:DemoFix3', {});
  }

  export async function readAsText(file: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ReadAsText', { file });
  }

  export async function readAsTextDapi(file: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ReadAsTextDapi', { file });
  }

  //Get biostructure by id as PDB
  export async function getBiostructureRcsbPdb(id: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:GetBiostructureRcsbPdb', { id });
  }

  //Get biostructure by id as mmCIF
  export async function getBiostructureRcsbMmcif(id: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:GetBiostructureRcsbMmcif', { id });
  }

  //Get biostructure by id as BinaryCIF
  export async function getBiostructureRcsbBcif(id: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:GetBiostructureRcsbBcif', { id });
  }

  //Packs BiostructureData value into JSON string
  export async function biostructureDataToJson(binary: boolean, data: any, ext: string, options: any): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:BiostructureDataToJson', { binary, data, ext, options });
  }

  export async function structure3D(molecule: any): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:Structure3D', { molecule });
  }

  //Extract protein sequences using Molstar parser
  export async function extractProteinSequencesMolstar(pdbId: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:ExtractProteinSequencesMolstar', { pdbId });
  }

  //Test the Molstar-based extraction
  export async function testMolstarExtraction(pdbId: string): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:TestMolstarExtraction', { pdbId });
  }

  //Test Molstar extraction with multiple PDBs
  export async function testMolstarMultiple(): Promise<any> {
    return await grok.functions.call('@datagrok/biostructure-viewer:TestMolstarMultiple', {});
  }
}
