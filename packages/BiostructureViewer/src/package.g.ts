import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: pdbCellRenderer
//output: grid_cell_renderer result
//meta.cellType: Molecule3D
//meta.columnTags: quality=Molecule3D
//meta.role: cellRenderer
export function Molecule3dCellRenderer() : any {
  return PackageFunctions.Molecule3dCellRenderer();
}

//name: chemCellRenderer
//output: grid_cell_renderer result
//meta.cellType: PDB_ID
//meta.role: cellRenderer
export function pdbIdCellRenderer() : any {
  return PackageFunctions.pdbIdCellRenderer();
}

//input: string pdbId 
export async function viewPdbById(pdbId: string) : Promise<void> {
  await PackageFunctions.viewPdbById(pdbId);
}

//input: string pdbData 
//input: string name 
export async function viewPdbByData(pdbData: string, name: string) : Promise<void> {
  await PackageFunctions.viewPdbByData(pdbData, name);
}

//output: object result
export function getNglGlService() : any {
  return PackageFunctions.getNglGlService();
}

//input: string content 
//input: string format 
//input: string name { optional: true }
export async function viewBiostructure(content: string, format?: string, name?: string) : Promise<void> {
  await PackageFunctions.viewBiostructure(content, format, name);
}

//description: Opens PDB file
//input: string fileContent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: mmcif, cifCore, pdb, gro
export async function importPdb(fileContent: string) : Promise<any> {
  return await PackageFunctions.importPdb(fileContent);
}

//description: Opens XYZ file
//input: string fileContent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: xyz
export async function importXYZ(fileContent: string) : Promise<any> {
  return await PackageFunctions.importXYZ(fileContent);
}

//description: Opens biostructure files supported with NGL
//input: string fileContent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: mmtf, cns, top, prmtop, ply, obj, ccp4
export async function importWithNgl(fileContent: string) : Promise<any> {
  return await PackageFunctions.importWithNgl(fileContent);
}

//description: Opens .pdbqt file with docking result ligand poses
//input: string fileContent 
//input: bool test = false { optional: true }
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: pdbqt
export async function importPdbqt(fileContent: string, test: boolean) : Promise<any> {
  return await PackageFunctions.importPdbqt(fileContent, test);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: mmtf,cns,top,prmtop,pqr
export function previewNglStructure(file: any) : any {
  return PackageFunctions.previewNglStructure(file);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: ply,obj
export function previewNglSurface(file: any) : any {
  return PackageFunctions.previewNglSurface(file);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: ccp4
export function previewNglDensity(file: any) : any {
  return PackageFunctions.previewNglDensity(file);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: mol,mol2,cif,mcif,mmcif,gro,pdb,pdbqt,ent,sd,xyz
export function previewBiostructureStructure(file: DG.FileInfo) : any {
  return PackageFunctions.previewBiostructureStructure(file);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: parm7,psf
export function previewBiostructureTopology(file: DG.FileInfo) : any {
  return PackageFunctions.previewBiostructureTopology(file);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: dsn6,brix,cube,cub,dx,dxbin,xplor,mrc,map
export function previewBiostructureDensity(file: DG.FileInfo) : any {
  return PackageFunctions.previewBiostructureDensity(file);
}

//input: file fi 
export async function openPdbResidues(fi: DG.FileInfo) : Promise<void> {
  await PackageFunctions.openPdbResidues(fi);
}

//name: PDB id viewer
//input: string pdbId { semType: PDB_ID }
//output: widget result
//meta.role: panel
export function pdbIdNglPanelWidget(pdbId: string) : any {
  return PackageFunctions.pdbIdNglPanelWidget(pdbId);
}

//description: Example app for NGL drawing in grid cells
export async function nglForGridTestApp() : Promise<void> {
  await PackageFunctions.nglForGridTestApp();
}

//description: Test app for NglViewer
export async function nglViewerApp() : Promise<void> {
  await PackageFunctions.nglViewerApp();
}

//description: Test app for BiostructureViewer (molstar)
export async function biostructureViewerApp() : Promise<void> {
  await PackageFunctions.biostructureViewerApp();
}

//description: Test app for BiotrackViewer (saguaro)
export async function biotrackViewerApp() : Promise<void> {
  await PackageFunctions.biotrackViewerApp();
}

//description: Test app for twin BiostructureViewer (molstar) and BiotrackViewer (saguaro)
export async function biostructureAndTrackViewerApp() : Promise<void> {
  await PackageFunctions.biostructureAndTrackViewerApp();
}

//name: ligandsWithNglApp
export async function ligandsWithNglApp() : Promise<void> {
  await PackageFunctions.ligandsWithNglApp();
}

//name: ligandsWithBiostructureApp
export async function ligandsWithBiostructureApp() : Promise<void> {
  await PackageFunctions.ligandsWithBiostructureApp();
}

//name: biostructureDataProviderApp
export async function biostructureDataProviderApp() : Promise<void> {
  await PackageFunctions.biostructureDataProviderApp();
}

//name: NGL
//description: 3D structure viewer for large biological molecules (proteins, DNA, and RNA)
//output: viewer result
//meta.keywords: PDB, Biostructure
//meta.icon: files/icons/ngl-viewer.svg
//meta.role: viewer,panel
export function nglViewer() {
  return PackageFunctions.nglViewer();
}

//name: Biostructure
//description: 3D structure molstar RCSB viewer for large biological molecules (proteins, DNA, and RNA)
//output: viewer result
//meta.keywords: Molstar, PDB
//meta.icon: files/icons/biostructure-viewer.svg
//meta.role: viewer,panel
export function molstarViewer() {
  return PackageFunctions.molstarViewer();
}

//name: Biotrack
//description: structure polymer annotation tracks
//output: viewer result
//meta.keywords: PDB, track
//meta.showInGallery: false
//meta.role: viewer,panel
export function saguaroViewer() {
  return PackageFunctions.saguaroViewer();
}

//output: object result
export async function getPdbHelper() : Promise<any> {
  return await PackageFunctions.getPdbHelper();
}

//name: dockingDemo
export async function dockingDemo() : Promise<void> {
  await PackageFunctions.dockingDemo();
}

//name: inGridDemo
export async function inGridDemo() : Promise<void> {
  await PackageFunctions.inGridDemo();
}

//name: Copy Biostructure raw value
//input: object gridCell 
export async function copyRawBiostructureValue(gridCell: any) : Promise<void> {
  await PackageFunctions.copyRawBiostructureValue(gridCell);
}

//name: Download Biostructure raw value
//input: object gridCell 
export async function downloadRawBiostructureValue(gridCell: any) : Promise<void> {
  await PackageFunctions.downloadRawBiostructureValue(gridCell);
}

//name: Show Biostructure Viewer menu item
//input: object gridCell 
export async function showBiostructureViewerMenuItem(gridCell: any) : Promise<void> {
  await PackageFunctions.showBiostructureViewerMenuItem(gridCell);
}

//name: Show NGL Viewer menu item
//input: object gridCell 
export async function showNglViewerMenuItem(gridCell: any) : Promise<void> {
  await PackageFunctions.showNglViewerMenuItem(gridCell);
}

//name: Open PDB residues table menu item
//input: object fi 
export async function openTableResiduesMenuItem(fi: DG.FileInfo) : Promise<void> {
  await PackageFunctions.openTableResiduesMenuItem(fi);
}

//description: Display ligand poses along the structure
//meta.demoPath: Bioinformatics | Docking Conformations
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Docking%20Conformations
//meta.demoWait: 3000
//test: demoBioDockingConformations() //wait: 3000 , timeout: 60000 
export async function demoBioDockingConformations() : Promise<void> {
  await PackageFunctions.demoBioDockingConformations();
}

//description: View structures PDB in grids
//meta.demoPath: Bioinformatics | Proteins
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Proteins
//meta.demoWait: 3000
//test: demoBioProteins() //wait: 3000 , timeout: 60000 
export async function demoBioProteins() : Promise<void> {
  await PackageFunctions.demoBioProteins();
}

//name: demoFix1
export async function demoFix1() : Promise<void> {
  await PackageFunctions.demoFix1();
}

//name: demoFix2
export async function demoFix2() : Promise<void> {
  await PackageFunctions.demoFix2();
}

//name: demoFix3
export async function demoFix3() : Promise<void> {
  await PackageFunctions.demoFix3();
}

//input: string file 
//output: string result
//meta.cache: client
//meta.cache.invalidateOn: 0 * * * *
export async function readAsText(file: string) : Promise<string> {
  return await PackageFunctions.readAsText(file);
}

//input: string file 
//output: string result
//meta.cache: client
//meta.cache.invalidateOn: 0 * * * *
export async function readAsTextDapi(file: string) : Promise<string> {
  return await PackageFunctions.readAsTextDapi(file);
}

//name: RCSB PDB
//description: Get biostructure by id as PDB
//input: string id 
//output: string result
//meta.dataProvider: Molecule3D
export async function getBiostructureRcsbPdb(id: string) : Promise<string> {
  return await PackageFunctions.getBiostructureRcsbPdb(id);
}

//name: RCSB mmCIF
//description: Get biostructure by id as mmCIF
//input: string id 
//output: string result
//meta.dataProvider: Molecule3D
//meta.cache: client
//meta.cache.invalidateOn: 0 * * * *
export async function getBiostructureRcsbMmcif(id: string) : Promise<string> {
  return await PackageFunctions.getBiostructureRcsbMmcif(id);
}

//name: RCSB bCIF
//description: Get biostructure by id as BinaryCIF
//input: string id = '1QBS' 
//output: string result
//meta.dataProvider: Molecule3D
//meta.cache: client
//meta.cache.invalidateOn: 0 * * * *
export async function getBiostructureRcsbBcif(id: string) : Promise<string> {
  return await PackageFunctions.getBiostructureRcsbBcif(id);
}

//description: Packs BiostructureData value into JSON string
//input: bool binary 
//input: object data 
//input: string ext 
//input: map options { optional: true }
//output: string result
export function biostructureDataToJson(binary: boolean, data: any, ext: string, options?: any) : string {
  return PackageFunctions.biostructureDataToJson(binary, data, ext, options);
}

//name: 3D Structure
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export function structure3D(molecule: DG.SemanticValue) : any {
  return PackageFunctions.structure3D(molecule);
}

//name: Fetch PDB Sequences
//description: For a user-selected table and PDB ID column, fetches protein sequences and adds them as new columns.
//input: dataframe table 
//input: column pdbColumn { semType: PDB_ID }
//top-menu: Bio | Transform | Fetch PDB Sequences...
export async function fetchSequencesFromPdb(table: DG.DataFrame, pdbColumn: DG.Column) : Promise<void> {
  await PackageFunctions.fetchSequencesFromPdb(table, pdbColumn);
}
