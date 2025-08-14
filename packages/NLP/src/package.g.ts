import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Translation
//input: file textfile 
//output: widget result
//condition: isTextFile(textfile)
export async function translationPanel(textfile: DG.FileInfo) {
  return PackageFunctions.translationPanel(textfile);
}

//name: exportFunc
//tags: init
export async function initAWS() {
  return PackageFunctions.initAWS();
}

//name: Compute Text Embeddings
//description: Compute text embeddings using UMAP
//top-menu: ML | Text Clustering...
export function computeEmbds() {
  return PackageFunctions.computeEmbds();
}

//name: Stem Column
//input: column col { semType: Text }
//input: string metric 
//input: int minimumCharactersCount { min: 0; max: 100; optional: true; default: 1 }
//output: object result
//meta.supportedSemTypes: Text
//meta.supportedDistanceFunctions: Common Items
export function stemColumnPreprocessingFunction(col: DG.Column, metric: string, minimumCharactersCount: number) {
  return PackageFunctions.stemColumnPreprocessingFunction(col, metric, minimumCharactersCount);
}

//name: Radial Coloring
//input: column col1 
//input: column col2 
export function radialColoring(col1: DG.Column, col2: DG.Column) {
  return PackageFunctions.radialColoring(col1, col2);
}

//name: Distance
//input: string query { semType: Text }
//output: widget result
//condition: true
export function distance(query: string) {
  return PackageFunctions.distance(query);
}

//name: Similar
//input: string query { semType: Text }
//output: widget result
//condition: true
export function similar(query: string) {
  return PackageFunctions.similar(query);
}

//name: Sentence Embeddings
//tags: dim-red-preprocessing-function
//input: column col { semType: Text }
//input: string metric 
//output: object result
//meta.supportedSemTypes: Text
//meta.supportedDistanceFunctions: Vector Cosine, Euclidean, Manhattan
export async function sentenceEmbeddingsPreprocessingFunction(col: DG.Column, metric: string) {
  return PackageFunctions.sentenceEmbeddingsPreprocessingFunction(col, metric);
}

//name: getEmbeddings
//input: list<string> sentences 
//output: string result
export async function getEmbeddings(sentences: string[]) {
  return PackageFunctions.getEmbeddings(sentences);
}

//name: Sentence Similarity Search
//tags: viewer
//output: viewer result
export function sentenceSearchViewer() {
  return PackageFunctions.sentenceSearchViewer();
}

//name: sentenceSearch
//output: viewer result
//top-menu: ML | Sentence Similarity Search
export function sentenceSearchTopMenu() {
  return PackageFunctions.sentenceSearchTopMenu();
}
