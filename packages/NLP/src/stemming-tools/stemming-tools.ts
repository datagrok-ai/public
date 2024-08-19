/* eslint-disable no-unused-vars */
/* eslint-disable valid-jsdoc */
// Stemming-based search (SBS) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {stemmer} from 'stemmer';
import {UMAP} from 'umap-js';
import {STR_METRIC_TYPE, MetricInfo, DISTANCE_TYPE, getSingleSourceMetricsMap, getDefaultDistnce,
  getDistanceFn, getMetricFn} from './metrics';
import {RATIO, INFTY, TINY, MIN_CHAR_COUNT, STOP_WORDS, SEPARATORS} from './constants';

/** Word information: index in the dictionary and frequency. */
type WordInfo = {
  index: number,
  frequency: number,
};

/** Cash of stemmed text column type. */
type StemCash = {
  dfName: string | undefined,
  colName: string | undefined,
  dictionary: Map<string, WordInfo> | undefined,
  indices: DG.Column | undefined,
  mostCommon: Set<number> | undefined,
  metricDef: Map<string, MetricInfo> | undefined,
  aggrDistance: DISTANCE_TYPE | undefined,
  filters: DG.FilterGroup | undefined,
};

/** Column (feature) usage specification: type, metric & use. */
export type ColUseInfo = {
  type: string,
  metric: MetricInfo,
  use: boolean,
};

/** UMAP Settings */
export type UmapSettings = {
  components: number,
  epochs: number,
  neighbors: number,
  minDist: number,
  spread: number,
};

/** Cash of stemmed text column. */
// eslint-disable-next-line no-var
export var stemCash: StemCash = {
  dfName: undefined,
  colName: undefined,
  dictionary: undefined,
  indices: undefined,
  mostCommon: undefined,
  metricDef: undefined,
  aggrDistance: undefined,
  filters: undefined,
};

/** Feature operating specification. */
type FeatureSpecification = {
  data: any,
  metricFn: (a: any, b: any) => number,
};

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
export function stemmColumn(source: DG.Column<string>, minCharCount: number):
{dict: Map<string, WordInfo>, indices: DG.Column, mostCommon: Set<number>} {
  const dict = new Map<string, WordInfo>();
  const size = source.length;
  const indices = DG.Column.fromType(DG.COLUMN_TYPE.OBJECT, `${source.name} (inds)`, size);
  const mostCommon = new Set<number>();

  for (let i = 0; i < size; ++i) {
    const words = getStemmedWords(source.get(i) ?? '', minCharCount);
    const inds = new Uint16Array(words.length);
    let j = 0;

    words.forEach((w) => {
      if (!dict.has(w)) {
        dict.set(w, {index: dict.size, frequency: 1});
        inds[j] = dict.size;
      } else {
        const info = dict.get(w);
        ++(info!.frequency);
        dict.set(w, info!);
        inds[j] = info!.index;
      }

      ++j;
    });

    indices.set(i, inds.sort());
  }

  const limit = RATIO * dict.size;

  for (const val of dict.values()) {
    if (val.frequency > limit)
      mostCommon.add(val.index);
  }

  return {dict: dict, indices: indices, mostCommon: mostCommon};
} // stemmColumn

/** Return number of common elements: the most common items are skipped. */
function commonItemsCount(arr1: Uint16Array, arr2: Uint16Array): number {
  const len1 = arr1.length;
  const len2 = arr2.length;
  let count = 0;
  let i1 = 0;
  let i2 = 0;

  while ((i1 < len1) && (i2 < len2)) {
    if (arr1[i1] === arr2[i2]) {
      if (!stemCash.mostCommon?.has(arr1[i1]))
        ++count;
      ++i1;
      ++i2;
    } else if (arr1[i1] < arr2[i2])
      ++i1;
    else
      ++i2;
  }

  return count;
}

/** TODO: to remove, old version*/
function distFn(arr1: number[], arr2: number[]): number {
  const idx1 = arr1[0];
  const idx2 = arr2[0];

  if (idx1 === idx2)
    return 0;

  const ref1 = stemCash.indices?.get(idx1) as Uint16Array;
  const ref2 = stemCash.indices?.get(idx2) as Uint16Array;

  const len1 = ref1.length;
  const len2 = ref2.length;

  if ((len1 === 0) || (len2 === 0))
    return INFTY;

  const common = commonItemsCount(ref1, ref2);

  if (common < 2)
  //if (common === 0)
    return INFTY;

  //return ((len1 < len2) ? len1 : len2) / common - 1;
  return ((len1 < len2) ? 1 : len2) / common - 1;
}

/** Return string with marked words and a dictionary of common words. */
export function getMarkedString(curIdx: number, queryIdx: number, strToBeMarked: string): HTMLElement[] {
  const df = grok.shell.t;
  const size = strToBeMarked.length;

  if (size === 0)
    return [];

  const marked = [] as any[];

  // stemmed words of the query text.
  const wordsOfQuery = stemCash.indices?.get(queryIdx);

  // indeces specifying the current subsequence
  let start = 0;
  let end = 0;

  // subsequence type: word or not
  let isWord = !SEPARATORS.includes(strToBeMarked[start]);

  const getWgt = (word: string, toMark: boolean) => {
    const p = ui.inlineText([word]);
    //@ts-ignore
    $(p).css({'cursor': 'pointer'});
    p.onclick = () => {
      df.currentCell = df.cell(curIdx, stemCash.colName!);
    };
    p.oncontextmenu = (e) => {
      e.stopImmediatePropagation();
      e.preventDefault();

      stemCash.filters = grok.shell.getTableView(df.name).getFiltersGroup();

      setTimeout(() => {
        const state = stemCash.filters!.getStates(stemCash.colName!, 'text')[0];
        // @ts-ignore
        state.gridNames.push(word);
          stemCash.filters!.updateOrAdd(state as DG.FilterState);
      }, 50);
    };
    ui.tooltip.bind(p, 'Click to navigate. Right-click to add to filters');

    if (toMark)
      //@ts-ignore
      $(p).css({'font-weight': 'bold'});

    return p;
  };

  // scan the string that should be marked
  while (start < size) {
    if (isWord) {
      while (!SEPARATORS.includes(strToBeMarked[end])) {
        ++end;

        if (end === size)
          break;
      }

      const word = strToBeMarked.slice(start, end);
      const wordLow = word.toLowerCase();

      if (!STOP_WORDS.includes(wordLow)) {
        const value = stemCash.dictionary?.get(stemmer(wordLow));

        if (wordsOfQuery.includes(value?.index))
          marked.push(getWgt(word, true));
        else
          marked.push(getWgt(word, false));
      } else
        marked.push(getWgt(word, false));
    // eslint-disable-next-line brace-style
    } // if (isWord)

    else {
      while (SEPARATORS.includes(strToBeMarked[end])) {
        ++end;

        if (end === size)
          break;
      }

      marked.push(getWgt(strToBeMarked.slice(start, end), false));
    }

    start = end;
    isWord = !isWord;
  } // while

  return marked;
} // getMarkedStringAndCommonWordsMap

/** Set stemming cash data. */
export function setStemmingCash(df: DG.DataFrame, source: DG.Column): void {
  if ((stemCash.dfName !== df.name) || (stemCash.colName !== source.name)) {
    stemCash.dfName = df.name;
    stemCash.colName = source.name;
    const stemmingRes = stemmColumn(source, MIN_CHAR_COUNT);
    stemCash.dictionary = stemmingRes.dict;
    stemCash.indices = stemmingRes.indices;
    stemCash.mostCommon = stemmingRes.mostCommon;
    stemCash.metricDef = getSingleSourceMetricsMap(source);
    stemCash.aggrDistance = getDefaultDistnce();
    stemCash.filters = undefined;
  }
}

/** Return text similarity distance function. */
function getTextDistFn(metric: MetricInfo): (query: Uint16Array, current: Uint16Array) => number {
  switch (metric.type) {
  case STR_METRIC_TYPE.STEMMING_BASED:
    return (query: Uint16Array, current: Uint16Array) => {
      const curLen = current.length;

      if (curLen === 0)
        return INFTY;

      return metric.weight * ((query.length < curLen) ? 1 : curLen) / (commonItemsCount(query, current) + TINY);
    };

  case STR_METRIC_TYPE.EQUALITY:
    return (query: Uint16Array, current: Uint16Array) => {
      if (query === current)
        return metric.weight;

      return 0;
    };

  default:
    throw new Error('Unsupported metric');
  }
}

/** Return closest elements indeces array. */
export function getClosest(df: DG.DataFrame, idx: number, count: number): Uint32Array {
  const indices = stemCash.indices!;
  const size = df.rowCount;
  const queryText = indices.get(idx) as Uint16Array;
  let textColDistFn = (a: Uint16Array, b: Uint16Array) => 0;
  let isTargetColUsed = false;
  const queryLen = queryText.length;

  // check sizes
  if ((size < count) || (queryLen === 0) || (stemCash.metricDef === undefined))
    return new Uint32Array([...Array(size).keys()]);

  // check metrics specification
  if (stemCash.metricDef.size < 1)
    return new Uint32Array([...Array(count).keys()]);

  const featuresCount = stemCash.metricDef.size;
  const featuresMetricVals = new Float32Array(featuresCount);
  const featureDistance = getDistanceFn(stemCash.aggrDistance!, featuresCount);
  const featureMap = new Map<string, FeatureSpecification>();
  const queryFeatures = new Map<string, any>();

  // define metrics
  for (const [key, value] of stemCash.metricDef) {
    if (key === stemCash.colName!) {
      isTargetColUsed = true;
      textColDistFn = getTextDistFn(value);
    } else {
      featureMap.set(key, {
        data: df.col(key)?.getRawData(),
        metricFn: getMetricFn(value),
      });

      queryFeatures.set(key, featureMap.get(key)?.data[idx]);
    }
  }

  const closestInds = new Uint32Array(count);
  const closestDists = new Float32Array(count);

  // basic initialization
  for (let i = 0; i < count; ++i) {
    let j = 0;

    for (const [key, value] of featureMap) {
      featuresMetricVals[j] = value.metricFn(queryFeatures.get(key), value.data[i]);
      ++j;
    }

    if (isTargetColUsed)
      featuresMetricVals[j] = textColDistFn(queryText, indices.get(i));

    closestInds[i] = i;
    closestDists[i] = featureDistance(featuresMetricVals);
  }

  let maxRelIdx = closestDists.reduce((r, v, i, a) => v <= a[r] ? r : i, -1);

  // search for the closest elements
  for (let i = count; i < size; ++i) {
    let j = 0;

    for (const [key, value] of featureMap) {
      featuresMetricVals[j] = value.metricFn(queryFeatures.get(key), value.data[i]);
      ++j;
    }

    if (isTargetColUsed)
      featuresMetricVals[j] = textColDistFn(queryText, indices.get(i));

    const curDist = featureDistance(featuresMetricVals);

    if (curDist < closestDists[maxRelIdx]) {
      closestDists[maxRelIdx] = curDist;
      closestInds[maxRelIdx] = i;
      maxRelIdx = closestDists.reduce((r, v, i, a) => v <= a[r] ? r : i, -1);
    }
  }

  return closestInds.filter((i) => i !== idx);
} // getClosest

/** Distance function fo the UMAP algorithm. */
function umapDistFn(arr1: number[], arr2: number[]): number {
  const df = grok.shell.table(stemCash.dfName!);

  const idx1 = arr1[0];
  const idx2 = arr2[0];

  const featuresCount = stemCash.metricDef!.size;
  const featuresMetricVals = new Float32Array(featuresCount);

  const i = 0;

  for (const [key, value] of stemCash.metricDef!) {
    if (key === stemCash.colName)
      featuresMetricVals[i] = getTextDistFn(value)(stemCash.indices!.get(idx1), stemCash.indices!.get(idx2));
    else
      featuresMetricVals[i] = getMetricFn(value)(df.get(key, idx1), df.get(key, idx2));
  }

  return getDistanceFn(stemCash.aggrDistance!, featuresCount)(featuresMetricVals);
}

/** Compute embeddings using the UMAP algorithm. */
export function getEmbeddings(table: DG.DataFrame, source: DG.Column, settings: UmapSettings): DG.Column[] {
  setStemmingCash(table, source);

  const data = [...Array(table.rowCount).keys()].map((i) => Array(10).fill(i));

  console.log(data);

  const umap = new UMAP({
    nComponents: settings.components,
    nEpochs: settings.epochs,
    nNeighbors: settings.neighbors,
    minDist: settings.minDist,
    spread: settings.spread,
    // @ts-ignore
    distanceFn: umapDistFn,
  });

  const embeddings = umap.fit(data!);

  const rowCount = embeddings.length;
  const range = [...Array(settings.components).keys()];
  const umapColumnsData = range.map((_) => new Float32Array(rowCount));

  for (let i = 0; i < rowCount; ++i) {
    for (let j = 0; j < settings.components; ++j)
      umapColumnsData[j][i] = embeddings[i][j];
  }

  return range.map((i) => DG.Column.fromFloat32Array('UMAP' + i.toString(), umapColumnsData[i]));
}
