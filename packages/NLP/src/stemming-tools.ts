// Stemming-based search (SBS) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {stemmer} from 'stemmer';
import { UMAP } from 'umap-js';

type StemBuffer = {
  dfName: string | undefined,
  colName: string | undefined,
  dictionary: Map<string, number> | undefined,
  indices: DG.Column | undefined,
};

export var stemBuffer: StemBuffer = {
  dfName: undefined,
  colName: undefined,
  dictionary: undefined,
  indices: undefined
};

const INFTY = 100000;
export const MIN_CHAR_COUNT = 1;

const STOP_WORDS = ["a", "about", "above", "after", "again", "against", "all", "am", "an", "and", "any", "are", "as", "at", 
  "be", "because", "been", "before", "being", "below", "between", "both", "but", "by", 
  "could", 
  "did", "do", "does", "doing", "down", "during", 
  "each", 
  "few", "for", "from", "further", 
  "had", "has", "have", "having", "he", "he'd", "he'll", "he's", "her", "here", "here's", "hers", "herself", "him", "himself", "his", "how", "how's", 
  "i", "i'd", "i'll", "i'm", "i've", "if", "in", "into", "is", "it", "it's", "its", "itself", "i.e",
  "let's", 
  "me", "more", "most", "my", "myself", 
  "nor", 'nan', 
  "of", "on", "once", "only", "or", "other", "ought", "our", "ours", "ourselves", "out", "over", "own", 
  "same", "she", "she'd", "she'll", "she's", "should", "so", "some", "such", 
  "than", "that", "that's", "the", "their", "theirs", "them", "themselves", "then", "there", "there's", "these", "they", "they'd", "they'll", "they're", "they've", "this", "those", "through", "to", "too", 
  "under", "until", "up", 
  "very", 
  "was", "we", "we'd", "we'll", "we're", "we've", "were", "what", "what's", "when", "when's", "where", "where's", "which", "while", "who", "who's", "whom", "why", "why's", "with", "would", 
  "you", "you'd", "you'll", "you're", "you've", "your", "yours", "yourself", "yourselves"];

const SEPARATORS = '[].,:;?!(){} \n\t#';

/** Get words of input string, stop words are emoved. */
function getWords(str: string, minCharCount: number): string[] {
  return str.toLowerCase()
    .replace(/[\s.,:;?!(){}#]/g, ' ')
    .split(' ')
    .filter((word) => ((word.length > minCharCount) && (word !== '\n') && (!STOP_WORDS.includes(word))));
}

/** Get an array of stemmed words. */
function getStemmedWords(str: string, minCharCount: number): string[] {
  return [...new Set(getWords(str, minCharCount).map((word) => stemmer(word)))];
}

/** Return dictionary of stemmed words & a column of indices. */
export function stemmColumn(source: DG.Column<string>, minCharCount: number): {dict: Map<string, number>, indices: DG.Column} {
  const dict = new Map<string, number>();
  const size = source.length;
  const indices = DG.Column.fromType(DG.COLUMN_TYPE.OBJECT, `${source.name} (inds)`, size);

  for (let i = 0; i < size; ++i) {
    const words = getStemmedWords(source.get(i) ?? '', minCharCount);
    const inds = new Uint16Array(words.length);
    let j = 0; 

    words.forEach((w) => {
      if (!dict.has(w))
        dict.set(w, dict.size);

      inds[j] = dict.get(w)!;
      
      ++j;
    });

    indices.set(i, inds.sort());
  }

  return {dict: dict, indices: indices};
}

/** Return number of common elements. */
function commonItemsCount(arr1: Uint16Array, arr2: Uint16Array): number {
  const len1 = arr1.length;
  const len2 = arr2.length;
  let count = 0;
  let i1 = 0;
  let i2 = 0;

  while((i1 < len1) && (i2 < len2))
    if (arr1[i1] === arr2[i2]) {
      ++count;
      ++i1;
      ++i2;
    }
    else if (arr1[i1] < arr2[i2])
      ++i1;
    else
      ++i2;
  
  return count;
}

/** Return indices of the closest items with respect to the overlap coefficient. */
export function getClosest(indices: DG.Column, idx: number, count: number): Uint32Array {
  const size = indices.length;
  const query = indices.get(idx) as Uint16Array;
  const queryLen = query.length;

  if ((size < count) || (queryLen === 0))
    return new Uint32Array([...Array(size).keys()]);

  const closestInds = new Uint32Array(count);
  const closestDists = new Float32Array(count);
    
  for (let i = 0; i < count; ++i) {
    const current = indices.get(i) as Uint16Array;
    const curLen = current.length;
 
    closestInds[i] = i;

    if (curLen === 0) 
      closestDists[i] = 0;    
    else 
      closestDists[i] = commonItemsCount(query, current) / ((queryLen < curLen) ? queryLen : curLen);
  }
    
  let minRelIdx = closestDists.reduce((r, v, i, a) => v >= a[r] ? r : i, -1);
  
  for (let i = count; i < size; ++i) {
    const current = indices.get(i) as Uint16Array;
    const curLen = current.length;

    if (curLen > 0) {      
      const curRelDist = commonItemsCount(query, current) / ((queryLen < curLen) ? queryLen : curLen);     

      if (curRelDist > closestDists[minRelIdx]) {
        closestDists[minRelIdx] = curRelDist;
        closestInds[minRelIdx] = i;
        minRelIdx = closestDists.reduce((r, v, i, a) => v >= a[r] ? r : i, -1);
      }
    }
  }
  
  return closestInds.filter(i => i !== idx);
}

function distaFn(arr1: Uint16Array, arr2: Uint16Array): number {
  const len1 = arr1.length;
  const len2 = arr2.length;

  if ((len1 === 0) || (len2 === 0))
    return INFTY;

  const common = commonItemsCount(arr1, arr2);

  if (common === 0)
    return INFTY;

  return ((len1 < len2) ? len1 : len2) / common;
}

function distFn(arr1: number[], arr2: number[]): number {
  const idx1 = arr1[0];
  const idx2 = arr2[0];

  if (idx1 === idx2)
    return 0;

  const ref1 = stemBuffer.indices?.get(idx1) as Uint16Array;
  const ref2 = stemBuffer.indices?.get(idx2) as Uint16Array;

  const len1 = ref1.length;
  const len2 = ref2.length;  

  if ((len1 === 0) || (len2 === 0))
    return INFTY;

  const common = commonItemsCount(ref1, ref2);  

  if (common < 2)
  //if (common === 0)
    return INFTY;

  return 1.0 / common - 1;
  //return ((len1 < len2) ? len1 : len2) / common - 1;
  //return ((len1 > len2) ? len1 : len2) / common - 1;
  //return (len1 + len2 - common) / common - 1;
}

/** Get embeddings. */
export function getEmbeddings(table: DG.DataFrame, source: DG.Column, components: number, epochs: number, neighbors: number, minDist: number, spread: number): DG.Column[] {
  if ((stemBuffer.dfName !== table.name) || (stemBuffer.colName !== source.name)) {
    console.log('Getting buffer...');

    stemBuffer.dfName = table.name;
    stemBuffer.colName = source.name;
    stemBuffer.indices = stemmColumn(source, 1).indices;
  }
  else
    console.log('We have already an appropriate buffer!');  

  const DIM = 10;

  const data = [...Array(stemBuffer.indices?.length).keys()].map(i => Array(DIM).fill(i));

  //console.log(data);

  const umap = new UMAP({
    nComponents: components,
    nEpochs: epochs,
    nNeighbors: neighbors,
    minDist: minDist,
    spread: spread,
// @ts-ignore
    distanceFn: distFn,
  });

  const embeddings = umap.fit(data!);

  const rowCount = embeddings.length;
  const range = [...Array(components).keys()];
  const umapColumnsData = range.map(_ => new Float32Array(rowCount));  

  for (let i = 0; i < rowCount; ++i)
    for (let j = 0; j < components; ++j)
      umapColumnsData[j][i] = embeddings[i][j];  

  return range.map(i => DG.Column.fromFloat32Array('UMAP' + i.toString(), umapColumnsData[i]));
}

/**  */
export function getMarkedString(queryIdx: number, strToBeMarked: string, minCharCount: number): any[] {
  const size = strToBeMarked.length;

  if (size === 0)
    return [];

  let result = [] as any[];

  const wordsOfQuery = stemBuffer.indices?.get(queryIdx);

  let start = 0;
  let end = 0;
  
  let isWord = !SEPARATORS.includes(strToBeMarked[start]);

  while (start < size) {


    if (isWord) {

      while (!SEPARATORS.includes(strToBeMarked[end])) {
        ++end;

        if (end === size)
          break;
      }
      
      const word = strToBeMarked.slice(start, end);

      if (!STOP_WORDS.includes(word)) {
        const stemmed = stemmer(word.toLowerCase());
        const value = stemBuffer.dictionary?.get(stemmed);

        if (wordsOfQuery.includes(value)) {
          let p = ui.inlineText([word]);
          //@ts-ignore
          $(p).css({"font-weight": "bold"});
          result.push(p);
        }
        else
          result.push(word);
      }
      else
      result.push(word);
    }

    else {
      while (SEPARATORS.includes(strToBeMarked[end])) {
        ++end;

        if (end === size)
          break;
      }

      result.push(strToBeMarked.slice(start, end));
    }   

    start = end;
    isWord = !isWord;
  } // while

  return result;
}
