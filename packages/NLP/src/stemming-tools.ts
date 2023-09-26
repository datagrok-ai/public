// Stemming-based search (SBS) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {stemmer} from 'stemmer';

type StemBuffer = {
  dfName: string | undefined,
  colName: string | undefined,
  indices: DG.Column | undefined,
};

export var stemBuffer: StemBuffer = {
  dfName: undefined,
  colName: undefined,
  indices: undefined
};

const STOP_WORDS = ["a", "about", "above", "after", "again", "against", "all", "am", "an", "and", "any", "are", "as", "at", 
  "be", "because", "been", "before", "being", "below", "between", "both", "but", "by", 
  "could", 
  "did", "do", "does", "doing", "down", "during", 
  "each", 
  "few", "for", "from", "further", 
  "had", "has", "have", "having", "he", "he'd", "he'll", "he's", "her", "here", "here's", "hers", "herself", "him", "himself", "his", "how", "how's", 
  "i", "i'd", "i'll", "i'm", "i've", "if", "in", "into", "is", "it", "it's", "its", "itself", 
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

/** Get words of input string, stop words are emoved. */
function getWords(str: string, minCharCount: number): string[] {
  return str.toLowerCase()
    .replace(/[\s.,%:?!(){}#]/g, ' ')
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
