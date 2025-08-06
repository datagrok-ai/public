import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: init
//tags: init
//output: dynamic result
export async function init() {
  return PackageFunctions.init();
}

//name: pdbCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: Molecule3D
//meta.columnTags: quality=Molecule3D
export function Molecule3dCellRenderer() {
  return PackageFunctions.Molecule3dCellRenderer();
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-PDB_ID
//output: grid_cell_renderer result
//meta.cellType: PDB_ID
export function pdbIdCellRenderer() {
  return PackageFunctions.pdbIdCellRenderer();
}

//name: viewPdbById
//input: string pdbId 
//output: dynamic result
export async function viewPdbById(pdbId: string) {
  return PackageFunctions.viewPdbById(pdbId);
}

//name: viewPdbByData
//input: string pdbData 
//input: string name 
export async function viewPdbByData(pdbData: string, name: string) {
  return PackageFunctions.viewPdbByData(pdbData, name);
}

//name: getNglGlService
//output: dynamic result
export function getNglGlService() {
  return PackageFunctions.getNglGlService();
}

//name: viewBiostructure
//input: string content 
//input: string format 
//input: string name { optional: true }
export async function viewBiostructure(content: string, format: string, name: string) {
  return PackageFunctions.viewBiostructure(content, format, name);
}

//name: importPdb
//description: Opens PDB file
//tags: file-handler
//input: string fileContent 
//output: list result
//meta.ext: mmcif, cifCore, pdb, gro
export async function importPdb(fileContent: string) {
  return PackageFunctions.importPdb(fileContent);
}

//name: importXYZ
//description: Opens XYZ file
//tags: file-handler
//input: string fileContent 
//output: list result
//meta.ext: xyz
export async function importXYZ(fileContent: string) {
  return PackageFunctions.importXYZ(fileContent);
}

//name: importWithNgl
//description: Opens biostructure files supported with NGL
//tags: file-handler
//input: string fileContent 
//output: list result
//meta.ext: mmtf, cns, top, prmtop, ply, obj, ccp4
export async function importWithNgl(fileContent: string) {
  return PackageFunctions.importWithNgl(fileContent);
}

//name: importPdbqt
//description: Opens .pdbqt file with docking result ligand poses
//tags: file-handler
//input: string fileContent 
//input: bool test { optional: true; default: false }
//output: list result
//meta.ext: pdbqt
export async function importPdbqt(fileContent: string, test: boolean) {
  return PackageFunctions.importPdbqt(fileContent, test);
}

//name: previewNglStructure
//tags: fileViewer
//input: dynamic file 
//output: view result
//meta.fileViewer: mmtf,cns,top,prmtop,pqr
export function previewNglStructure(file: any) {
  return PackageFunctions.previewNglStructure(file);
}

//name: previewNglSurface
//tags: fileViewer
//input: dynamic file 
//output: dynamic result
//meta.fileViewer: ply,obj
export function previewNglSurface(file: any) {
  return PackageFunctions.previewNglSurface(file);
}

//name: previewNglDensity
//tags: fileViewer
//input: dynamic file 
//output: dynamic result
//meta.fileViewer: ccp4
export function previewNglDensity(file: any) {
  return PackageFunctions.previewNglDensity(file);
}

//name: previewBiostructureStructure
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: mol,mol2,cif,mcif,mmcif,gro,pdb,pdbqt,ent,sd,xyz
export function previewBiostructureStructure(file: DG.FileInfo) {
  return PackageFunctions.previewBiostructureStructure(file);
}

//name: previewBiostructureTopology
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: parm7,psf
export function previewBiostructureTopology(file: DG.FileInfo) {
  return PackageFunctions.previewBiostructureTopology(file);
}

//name: previewBiostructureDensity
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: dsn6,brix,cube,cub,dx,dxbin,xplor,mrc,map
export function previewBiostructureDensity(file: DG.FileInfo) {
  return PackageFunctions.previewBiostructureDensity(file);
}

//name: openPdbResidues
//input: file fi 
export async function openPdbResidues(fi: DG.FileInfo) {
  return PackageFunctions.openPdbResidues(fi);
}

//name: PDB id viewer
//tags: panel
//input: string pdbId { semType: PDB_ID }
//output: widget result
export function pdbIdNglPanelWidget(pdbId: string) {
  return PackageFunctions.pdbIdNglPanelWidget(pdbId);
}

//name: nglForGridTestApp
//description: Example app for NGL drawing in grid cells
//output: dynamic result
export async function nglForGridTestApp() {
  return PackageFunctions.nglForGridTestApp();
}

//name: nglViewerApp
//description: Test app for NglViewer
//output: dynamic result
export async function nglViewerApp() {
  return PackageFunctions.nglViewerApp();
}

//name: biostructureViewerApp
//description: Test app for BiostructureViewer (molstar)
export async function biostructureViewerApp() {
  return PackageFunctions.biostructureViewerApp();
}

//name: biotrackViewerApp
//description: Test app for BiotrackViewer (saguaro)
export async function biotrackViewerApp() {
  return PackageFunctions.biotrackViewerApp();
}

//name: biostructureAndTrackViewerApp
//description: Test app for twin BiostructureViewer (molstar) and BiotrackViewer (saguaro)
export async function biostructureAndTrackViewerApp() {
  return PackageFunctions.biostructureAndTrackViewerApp();
}

//name: ligandsWithNglApp
export async function ligandsWithNglApp() {
  return PackageFunctions.ligandsWithNglApp();
}

//name: ligandsWithBiostructureApp
export async function ligandsWithBiostructureApp() {
  return PackageFunctions.ligandsWithBiostructureApp();
}

//name: biostructureDataProviderApp
export async function biostructureDataProviderApp() {
  return PackageFunctions.biostructureDataProviderApp();
}

//name: NGL
//description: 3D structure viewer for large biological molecules (proteins, DNA, and RNA)
//tags: panel, viewer
//output: viewer result
//meta.keywords: PDB, Biostructure
//meta.icon: files/icons/ngl-viewer.svg
export function nglViewer() {
  return PackageFunctions.nglViewer();
}

//name: Biostructure
//description: 3D structure molstar RCSB viewer for large biological molecules (proteins, DNA, and RNA)
//tags: panel, viewer
//output: viewer result
//meta.keywords: Molstar, PDB
//meta.icon: files/icons/biostructure-viewer.svg
export function molstarViewer() {
  return PackageFunctions.molstarViewer();
}

//name: Biotrack
//description: structure polymer annotation tracks
//tags: panel, viewer
//output: viewer result
//meta.keywords: PDB, track
export function saguaroViewer() {
  return PackageFunctions.saguaroViewer();
}

//name: getPdbHelper
//output: dynamic result
export async function getPdbHelper() {
  return PackageFunctions.getPdbHelper();
}

//name: dockingDemo
//output: dynamic result
export async function dockingDemo() {
  return PackageFunctions.dockingDemo();
}

//name: inGridDemo
//output: dynamic result
export async function inGridDemo() {
  return PackageFunctions.inGridDemo();
}

//name: Copy Biostructure raw value
//input: dynamic gridCell 
export async function copyRawBiostructureValue(gridCell: any) {
  return PackageFunctions.copyRawBiostructureValue(gridCell);
}

//name: Download Biostructure raw value
//input: dynamic gridCell 
export async function downloadRawBiostructureValue(gridCell: any) {
  return PackageFunctions.downloadRawBiostructureValue(gridCell);
}

//name: Show Biostructure Viewer menu item
//input: dynamic gridCell 
//output: dynamic result
export async function showBiostructureViewerMenuItem(gridCell: any) {
  return PackageFunctions.showBiostructureViewerMenuItem(gridCell);
}

//name: Show NGL Viewer menu item
//input: dynamic gridCell 
//output: dynamic result
export async function showNglViewerMenuItem(gridCell: any) {
  return PackageFunctions.showNglViewerMenuItem(gridCell);
}

//name: Open PDB residues table menu item
//input: file fi 
//output: dynamic result
export async function openTableResiduesMenuItem(fi: DG.FileInfo) {
  return PackageFunctions.openTableResiduesMenuItem(fi);
}

//name: demoBioDockingConformations
//description: Display ligand poses along the structure
//meta.demoPath: Bioinformatics | Docking Conformations
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Docking%20Conformations
//meta.demoWait: 3000
//meta.demoSkip: GROK-15250
//test: demoBioDockingConformations() //wait: 3000 , timeout: 60000 , skip: GROK-15250 
export async function demoBioDockingConformations() {
  return PackageFunctions.demoBioDockingConformations();
}

//name: demoBioProteins
//description: View structures PDB in grids
//meta.demoPath: Bioinformatics | Proteins
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Proteins
//meta.demoWait: 3000
//meta.demoSkip: GROK-15250
//test: demoBioProteins() //wait: 3000 , timeout: 60000 , skip: GROK-15250 
export async function demoBioProteins() {
  return PackageFunctions.demoBioProteins();
}

//name: demoFix1
export async function demoFix1() {
  return PackageFunctions.demoFix1();
}

//name: demoFix2
export async function demoFix2() {
  return PackageFunctions.demoFix2();
}

//name: demoFix3
export async function demoFix3() {
  return PackageFunctions.demoFix3();
}

//name: readAsText
//input: string file 
//output: string result
//meta.cache: client
//meta.cacheInvalidateOn: 0 * * * *
export async function readAsText(file: string) {
  return PackageFunctions.readAsText(file);
}

//name: readAsTextDapi
//input: string file 
//output: dynamic result
//meta.cache: client
//meta.cacheInvalidateOn: 0 * * * *
export async function readAsTextDapi(file: string) {
  return PackageFunctions.readAsTextDapi(file);
}

//name: RCSB PDB
//description: Get biostructure by id as PDB
//input: string id 
//output: string result
//meta.dataProvider: Molecule3D
export async function getBiostructureRcsbPdb(id: string) {
  return PackageFunctions.getBiostructureRcsbPdb(id);
}

//name: RCSB mmCIF
//description: Get biostructure by id as mmCIF
//input: string id 
//output: string result
//meta.dataProvider: Molecule3D
//meta.cache: client
//meta.cacheInvalidateOn: 0 * * * *
export async function getBiostructureRcsbMmcif(id: string) {
  return PackageFunctions.getBiostructureRcsbMmcif(id);
}

//name: RCSB bCIF
//description: Get biostructure by id as BinaryCIF
//input: string id { default: 1QBS }
//output: string result
//meta.dataProvider: Molecule3D
//meta.cache: client
//meta.cacheInvalidateOn: 0 * * * *
export async function getBiostructureRcsbBcif(id: string) {
  return PackageFunctions.getBiostructureRcsbBcif(id);
}

//name: biostructureDataToJson
//description: Packs BiostructureData value into JSON string
//input: bool binary 
//input: dynamic data 
//input: string ext 
//input: dynamic options { optional: true }
//output: string result
export function biostructureDataToJson(binary: boolean, data: any, ext: string, options: any) {
  return PackageFunctions.biostructureDataToJson(binary, data, ext, options);
}

//name: 3D Structure
//tags: panel, bio, widgets
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
export function structure3D(molecule: DG.SemanticValue) {
  return PackageFunctions.structure3D(molecule);
}

//name: Fetch PDB Sequences
//description: For a user-selected table and PDB ID column, fetches protein sequences and adds them as new columns.
//input: dataframe table 
//input: column pdbColumn { semType: PDB_ID }
//top-menu: Bio | Transform | Fetch PDB Sequences...
export async function fetchSequencesFromPdb(table: DG.DataFrame, pdbColumn: DG.Column) {
  return PackageFunctions.fetchSequencesFromPdb(table, pdbColumn);
}
