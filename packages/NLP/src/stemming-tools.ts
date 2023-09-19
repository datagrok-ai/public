// Stemming-based search (SBS) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {stemmer} from 'stemmer';

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

export const STEMMED_SUFIX = ' (stemmed)';
const INDEX_COL_NAME = '#';

type Metrics = {
  absolute: number,
  relative: number
};

type Element = {
  idx: number,
  metrics: Metrics
};

type NumArrElement = {
  idx: number,
  val: number
};

type Closest = {
  abs: Element[],
  rel: Element[],
};

export type CorrectClassesCount = {
  absolute: number,
  absoluteBest: boolean,
  relative: number,
  relativeBest: boolean,
};

export type SearchResults = {
  indeces: number[],
  strings: string[]
};

export type SearchResultsDFs = {
  absClosest: DG.DataFrame, 
  relClosest: DG.DataFrame
};

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

/** Return common elements of string arrays. */
function commonElements(words1: string[], words2: string[]): string[] {  
  //return words1.filter((word) => words2.includes(word));  
  //return (words1.length < words2.length) ? words1.filter((word) => words2.includes(word)) : words2.filter((word) => words1.includes(word));
  return (words1.length > words2.length) ? words1.filter((word) => words2.includes(word)) : words2.filter((word) => words1.includes(word));
}
// Q.csv, 20117 rows:
/* Basic 
     before:        
         idx 587:  674.5 ms.
           idx 2:  675.3 ms.
       idx 13367:  878.85 ms. 
       
     "in min-length search":
         idx 587:  673.65 ms.
           idx 2:  681 ms.
       idx 13367:  879.25 ms.

     "in max-length search": !!!!
         idx 587:  674.3 ms.
           idx 2:  678.7 ms.
       idx 13367:  795.1 ms.
*/
/* Stemmed
     before:        
         idx 587:  31.2 ms.
           idx 2:  32.45 ms.
       idx 13367:  219.9 ms.
       
     "in min-length search":
         idx 587:  33.55 ms.
           idx 2:  35.8 ms.
       idx 13367:  222.55 ms.

     "in max-length search": !!!!
         idx 587:  31.7 ms.
           idx 2:  34 ms.
       idx 13367:  135.2 ms.
*/

/** Return common stemmed words of two strings. */
function commonStemmedWords(str1: string, str2: string, minCharCount: number): string[] {
  return commonElements(getStemmedWords(str1, minCharCount), getStemmedWords(str2, minCharCount));
}

/** Return absolute & relative stemming-based metrics. */
function distance(words1: string[], words2: string[]): Metrics {
  const divisor = words2.length;
  const absMetric = commonElements(words1, words2).length;
  
  return {
    absolute: absMetric, 
    relative: (divisor > 0) ? absMetric / divisor: 0
  }
}

/** Return array of metrics between the query and each element of the column base. Use only when processing small columns. */
function computeMetrics(query: string, base: DG.Column<string>, minCharCount: number): Element[] {
  const result = [] as Element[];
  const size = base.length;
  const stemmedQuery = getStemmedWords(query, minCharCount);

  for (let i = 0; i < size; ++i)
    result.push({
      idx: i,
      metrics: distance(
        stemmedQuery,
        getStemmedWords(base.get(i) ?? '', minCharCount)
    )});

  return result;
}

/** Return dataframe with stemming-based search result. */
function searchResultsDF(query: string, baseColumn: DG.Column<string>, closest: Element[]): DG.DataFrame {
  const closestItems = [];
  const closestNo = [];
  const closestAbsMetrics = [];
  const closestRelMetrics = [];
  const commonWords = [];
  
  for (const item of closest) {
    const current = baseColumn.get(item.idx) ?? '';
    closestNo.push(item.idx);
    closestItems.push(current);
    closestAbsMetrics.push(item.metrics.absolute);
    closestRelMetrics.push(item.metrics.relative);
    commonWords.push(commonStemmedWords(query, current, 1).join(', '));
  }
  
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.INT, INDEX_COL_NAME, closestNo),
    DG.Column.fromStrings(baseColumn.name, closestItems),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Metrics (absolute)', closestAbsMetrics),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Metrics (relative)', closestRelMetrics),
    DG.Column.fromStrings('Common words (stemmed)', commonWords)
  ]);
}

/** Return element of the array that has the smallest field specified by getField-function. */
function getMin(arr: Element[], getFieldFn: (elem: Element) => number): NumArrElement {
  const res = {idx: 0, val: getFieldFn(arr[0])};
  const size = arr.length;
  
  for (let i = 1; i < size; ++i) {
    const cur = getFieldFn(arr[i]);

    if (cur < res.val) {
      res.idx = i;
      res.val = cur;
    }
  }

  return res;
}

/** Return relative metric field. */
function getRel(elem: Element): number { return elem.metrics.relative; }

/** Return absolute metric field. */
function getAbs(elem: Element): number { return elem.metrics.absolute; }

/** Return two arrays of elements closest with respect to absolute & relative metrics. */
function closestElements(query: string, base: DG.Column<string>, closestCount: number, minCharCount: number): Closest {
  const size = base.length;

  if (size <= closestCount) {
    const metrics = computeMetrics(query, base, minCharCount);

    return {
      abs: [...metrics.sort((a, b) => b.metrics.absolute - a.metrics.absolute)],
      rel: [...metrics.sort((a, b) => b.metrics.relative - a.metrics.relative)]
    };
  }

  const closestAbs = [] as Element[];
  const closestRel = [] as Element[];
  const stemmedQuery = getStemmedWords(query, minCharCount);

  for (let i = 0; i < closestCount; ++i) {
    const elem = {idx: i, metrics: distance(stemmedQuery, getStemmedWords(base.get(i) ?? '', minCharCount))};
    closestAbs.push(elem);
    closestRel.push(elem);
  }

  let minAbs = getMin(closestAbs, getAbs);
  let minRel = getMin(closestRel, getRel);

  for (let i = closestCount; i < size; ++i) {
    const elem = {idx: i, metrics: distance(stemmedQuery, getStemmedWords(base.get(i) ?? '', minCharCount))};

    if (elem.metrics.absolute > minAbs.val) {
      closestAbs[minAbs.idx] = elem;
      minAbs = getMin(closestAbs, getAbs);
    }

    if (elem.metrics.relative > minRel.val) {
      closestRel[minRel.idx] = elem;
      minRel = getMin(closestRel, getRel);
    }
  }

  return {
    abs: [...closestAbs.sort((a, b) => b.metrics.absolute - a.metrics.absolute)],
    rel: [...closestRel.sort((a, b) => b.metrics.relative - a.metrics.relative)]
  };
}

/** Return two dataframes with stemming-based search results, obtained with respect to absolute & relative metrics. */
export function closestElementsDF(query: string, base: DG.Column<string>, closestCount: number, minCharCount: number): SearchResultsDFs {
  const closest = closestElements(query, base, closestCount, minCharCount);

  return {
    absClosest: searchResultsDF(query, base, closest.abs),
    relClosest: searchResultsDF(query, base, closest.rel)
  };
}

/** Return stemming-based search results: indices & strings. */
export function getSearchResults(query: string, base: DG.Column<string>, closestCount: number, minCharCount: number): SearchResults {
  const closest = closestElements(query, base, closestCount, minCharCount); 

  let indeces = closest.rel.map((item) => item.idx);
  indeces = indeces.concat(closest.abs.map((item) => item.idx).filter((idx) => !indeces.includes(idx)));
  
  return { 
    indeces: indeces, 
    strings: indeces.map((idx) => base.get(idx) ?? '') 
  };
}

/** Return column of stemmed words. */
export function stemmedColumn(col: DG.Column<string>, minCharCount: number): DG.Column<string> {
  const size = col.length;
  const stemmed = DG.Column.fromType(DG.COLUMN_TYPE.STRING, `${col.name}${STEMMED_SUFIX}`, size);

  for (let i = 0; i < size; ++i)
    stemmed.set(
      i,
      getStemmedWords(col.get(i) ?? '', minCharCount).join(' ')
    );

  return stemmed;
}

// HEURISTIC VALIDATION TOOLS

/** Return categories corresponding to the closest items */
export function getCategories(closest: DG.DataFrame, categories: DG.Column): DG.Column {
  const size = closest.rowCount;
  const res = DG.Column.fromType(DG.COLUMN_TYPE.STRING, categories.name, size);
  const indices = closest.col(INDEX_COL_NAME);

  for (let i = 0; i < size; ++i)
    res.set(i, categories.get(indices?.get(i)));

  return res; 
}

/** Find the number of occurrences of an element in a column. */
function countOccurs(element: any, col: DG.Column): number {  
  let res = 0;
  const size = col.length;

  for (let i = 0; i < size; ++i)
    if (col.get(i) === element)
      ++res;

  return res;
}

/** Heuristic validation. */
export function getValidationResults(idx: number, strings: DG.Column, categories: DG.Column, closestCount: number): CorrectClassesCount {
  const query = strings.get(((idx >= 0 ) && (idx < strings.length))? idx : 0) ?? '';
  const searchResults = closestElementsDF(query, strings, closestCount, 1);

  const absClasses = getCategories(searchResults.absClosest, categories);
  const relClasses = getCategories(searchResults.relClosest, categories);

  searchResults.absClosest.columns.add(absClasses);
  searchResults.relClosest.columns.add(relClasses);

  const queryClass = categories.get(idx);

  return {
    absolute: countOccurs(queryClass, absClasses) - 1,
    absoluteBest: queryClass === absClasses.get(1),
    relative: countOccurs(queryClass, relClasses) - 1,
    relativeBest: queryClass === relClasses.get(1),
  };
}

// EXPERIMENTAL TOOLS: discussion & optimization are required!

function computeMetricsUsingStemmedData(query: string, stemmed: DG.Column<string>, minCharCount: number): Element[] {
  const result = [] as Element[];
  const size = stemmed.length;
  const stemmedQuery = getStemmedWords(query, minCharCount);

  for (let i = 0; i < size; ++i)
    result.push({
      idx: i,
      metrics: distance(
        stemmedQuery,
        (stemmed.get(i) ?? '').split(' ')
    )});

  return result;
}

export function getClosestElementsDFusingStemmedData(query: string, stemmed: DG.Column<string>, closestCount: number, minCharCount: number): SearchResultsDFs {  
  const metrics = computeMetricsUsingStemmedData(query, stemmed, minCharCount);
  
  const closestAbs = [...metrics.sort((a, b) => b.metrics.absolute - a.metrics.absolute).slice(0, closestCount)];
  const closestRel = [...metrics.sort((a, b) => b.metrics.relative - a.metrics.relative).slice(0, closestCount)];
    
  return {
    absClosest: searchResultsDF(query, stemmed, closestAbs), 
    relClosest: searchResultsDF(query, stemmed, closestRel)
  };  
}

export function getClosestUsingStemmedData(query: string, source: DG.Column<string>, stemmed: DG.Column<string>, closestCount: number, minCharCount: number): SearchResults {
  const metrics = computeMetricsUsingStemmedData(query, stemmed, minCharCount);
    
  const closestRel = [...metrics.sort((a, b) => b.metrics.relative - a.metrics.relative).slice(0, closestCount)];
  const closestAbs = [...metrics.sort((a, b) => b.metrics.absolute - a.metrics.absolute).slice(0, 2 * closestCount)];

  let closestIndeces = closestRel.map((item) => item.idx);
  closestIndeces = closestIndeces.concat(closestAbs.map((item) => item.idx).filter((idx) => !closestIndeces.includes(idx)).slice(0, closestCount - 1));
  
  return {
    indeces: closestIndeces,
    strings: closestIndeces.map((idx) => source.get(idx) ?? '') 
  };
}
