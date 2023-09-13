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

type Metrics = {
  absolute: number,
  relative: number
};

type Element = {
  idx: number,
  metrics: Metrics
};

export type CorrectClassesCount = {
  absolute: number,
  absoluteBest: boolean,
  relative: number,
  relativeBest: boolean,
};

export type SearchResults = {
  absClosestIndeces: number[],
  absClosestStrings: string[],
  relClosestIndeces: number[],
  relClosestStrings: string[]
};

export type SearchResultsSplitted = {
  closestIndeces: number[],
  closestStrings: string[]
};

export type FullSearchResults = {
  absClosest: DG.DataFrame, 
  relClosest: DG.DataFrame
};

function getWords(str: string, minCharCount: number): string[] {
  return str.toLowerCase()
    .replace(/[\s.,%:?!(){}#]/g, ' ')
    .split(' ')
    .filter((word) => ((word.length > minCharCount) && (word !== '\n') && (!STOP_WORDS.includes(word))));
}

function getStemmedWords(str: string, minCharCount: number): string[] {
  return [...new Set(getWords(str, minCharCount).map((word) => stemmer(word)))];
}

function commonElements(words1: string[], words2: string[]): string[] {
  return words1.filter((word) => words2.includes(word));
}

function distance(words1: string[], words2: string[]): Metrics {
  const divisor = words2.length;
  const absMetric = commonElements(words1, words2).length;
  
  return {
    absolute: absMetric, 
    relative: (divisor > 0) ? absMetric / divisor: 0
  }
}

function commonStemmedWords(str1: string, str2: string, minCharCount: number): string[] {
  return commonElements(getStemmedWords(str1, minCharCount), getStemmedWords(str2, minCharCount));
}

function getMetrics(query: string, base: string[], minCharCount: number): Element[] {
  const result = [] as Element[];
  const size = base.length;

  for (let i = 0; i < size; ++i)
    result.push({
      idx: i,
      metrics: distance(
        getStemmedWords(query, minCharCount),
        getStemmedWords(base[i] ?? '', minCharCount)
    )});

  return result;
}

function computeMetrics(query: string, base: DG.Column<string>, minCharCount: number): Element[] {
  const result = [] as Element[];
  const size = base.length;

  for (let i = 0; i < size; ++i)
    result.push({
      idx: i,
      metrics: distance(
        getStemmedWords(query, minCharCount),
        getStemmedWords(base.get(i) ?? '', minCharCount)
    )});

  return result;
}

function getSearchResultsDataframe(query: string, baseColumn: DG.Column<string>, closest: Element[]): DG.DataFrame {
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
    DG.Column.fromList(DG.COLUMN_TYPE.INT, '#', closestNo),
    DG.Column.fromStrings(baseColumn.name, closestItems),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Metrics (absolute)', closestAbsMetrics),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Metrics (relative)', closestRelMetrics),
    DG.Column.fromStrings('Common words (stemmed)', commonWords)
  ]);
}

export function getClasses(results: DG.DataFrame, section: DG.Column): DG.Column {
  return DG.Column.fromList(
    DG.COLUMN_TYPE.STRING,
    section.name,
    results.columns.toList()[0].toList().map((idx) => section.get(idx))
  );
}

function count(element: any, col: DG.Column): number { return col.toList().filter((item) => (item === element)).length; }

export function getValidationResults(idx: number, questions: DG.Column, section: DG.Column, closestCount: number): CorrectClassesCount {
  const query = questions.get(((idx >= 0 ) && (idx < questions.length))? idx : 0) ?? '';
  const searchResults = getFullSearchResults(query, questions, closestCount, 1);

  const absClasses = getClasses(searchResults.absClosest, section);
  const relClasses = getClasses(searchResults.relClosest, section);

  searchResults.absClosest.columns.add(absClasses);
  searchResults.relClosest.columns.add(relClasses);

  const queryClass = section.get(idx);

  return {
    absolute: count(queryClass, absClasses) - 1,
    absoluteBest: queryClass === absClasses.get(1),
    relative: count(queryClass, relClasses) - 1,
    relativeBest: queryClass === relClasses.get(1),
  };
}

function getSearchResults(query: string, baseColumn: DG.Column<string>, closestCount: number, minCharCount: number): SearchResults {
  const base = baseColumn.toList() as string[];
  const metrics = getMetrics(query, base, minCharCount);
  
  const closestAbs = [...metrics.sort((a, b) => b.metrics.absolute - a.metrics.absolute).slice(0, closestCount)];
  const closestRel = [...metrics.sort((a, b) => b.metrics.relative - a.metrics.relative).slice(0, closestCount)];
  
  return {
    absClosestIndeces: closestAbs.map((item) => item.idx),
    absClosestStrings: closestAbs.map((item) => item.idx).map((idx) => base[idx]),
    relClosestIndeces: closestRel.map((item) => item.idx), 
    relClosestStrings: closestRel.map((item) => item.idx).map((idx) => base[idx])
  };
}

export function getSearchResultsSplitted(query: string, baseColumn: DG.Column<string>, closestCount: number, minCharCount: number): SearchResultsSplitted {
  const metrics = computeMetrics(query, baseColumn, minCharCount);
    
  const closestRel = [...metrics.sort((a, b) => b.metrics.relative - a.metrics.relative).slice(0, closestCount)];
  const closestAbs = [...metrics.sort((a, b) => b.metrics.absolute - a.metrics.absolute).slice(0, 2 * closestCount)];

  let closestIndeces = closestRel.map((item) => item.idx);
  closestIndeces = closestIndeces.concat(closestAbs.map((item) => item.idx).filter((idx) => !closestIndeces.includes(idx)).slice(0, closestCount - 1));
  
  return {
    closestIndeces: closestIndeces,
    closestStrings: closestIndeces.map((idx) => baseColumn.get(idx) ?? '') 
  };
}

export function getFullSearchResults(query: string, baseColumn: DG.Column<string>, closestCount: number, minCharCount: number): FullSearchResults {  
  const metrics = computeMetrics(query, baseColumn, minCharCount);
  
  const closestAbs = [...metrics.sort((a, b) => b.metrics.absolute - a.metrics.absolute).slice(0, closestCount)];
  const closestRel = [...metrics.sort((a, b) => b.metrics.relative - a.metrics.relative).slice(0, closestCount)];
    
  return {
    absClosest: getSearchResultsDataframe(query, baseColumn, closestAbs), 
    relClosest: getSearchResultsDataframe(query, baseColumn, closestRel)
  };  
}