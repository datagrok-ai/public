import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Translation
//input: file textfile 
//output: widget result
//condition: isTextFile(textfile)
export async function translationPanel(textfile: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.translationPanel(textfile);
}

//name: exportFunc
//meta.role: init
export async function initAWS() : Promise<void> {
  await PackageFunctions.initAWS();
}

//name: Compute Text Embeddings
//description: Compute text embeddings using UMAP
//top-menu: ML | Text Clustering...
export function computeEmbds() : void {
  PackageFunctions.computeEmbds();
}

//name: Stem Column
//input: column col { semType: Text }
//input: string metric 
//input: int minimumCharactersCount = 1 { min: 0; max: 100; optional: true }
//output: object result
//meta.supportedSemTypes: Text
//meta.supportedDistanceFunctions: Common Items
export function stemColumnPreprocessingFunction(col: DG.Column, metric: string, minimumCharactersCount: number) {
  return PackageFunctions.stemColumnPreprocessingFunction(col, metric, minimumCharactersCount);
}

//name: Radial Coloring
//input: column col1 
//input: column col2 
export function radialColoring(col1: DG.Column, col2: DG.Column) : void {
  PackageFunctions.radialColoring(col1, col2);
}

//name: Distance
//input: string query { semType: Text }
//output: widget result
//condition: true
export function distance(query: string) : any {
  return PackageFunctions.distance(query);
}

//name: Similar
//input: string query { semType: Text }
//output: widget result
//condition: true
export function similar(query: string) : any {
  return PackageFunctions.similar(query);
}

//name: Sentence Embeddings
//input: column col { semType: Text }
//input: string metric 
//output: object result
//meta.supportedSemTypes: Text
//meta.supportedDistanceFunctions: Vector Cosine, Euclidean, Manhattan
//meta.role: dimRedPreprocessingFunction
export async function sentenceEmbeddingsPreprocessingFunction(col: DG.Column, metric: string) {
  return await PackageFunctions.sentenceEmbeddingsPreprocessingFunction(col, metric);
}

//input: list<string> sentences 
//output: string result
export async function getEmbeddings(sentences: string[]) : Promise<any> {
  return await PackageFunctions.getEmbeddings(sentences);
}

//name: Sentence Similarity Search
//output: viewer result
//meta.showInGallery: false
//meta.role: viewer
export function sentenceSearchViewer() : any {
  return PackageFunctions.sentenceSearchViewer();
}

//name: sentenceSearch
//output: viewer result
//top-menu: ML | Sentence Similarity Search
export function sentenceSearchTopMenu() {
  return PackageFunctions.sentenceSearchTopMenu();
}
