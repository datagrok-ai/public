// Stemming-based search (SBS) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {stemmer} from 'stemmer';
import { UMAP } from 'umap-js';
import {STR_METRIC_TYPE, MetricType, MetricInfo, getDefaultMetric, 
  getChoicesList, DISTANCE_TYPE, getOneHotMetricsMap, getDefaultDistnce, DISTANCE_TYPES,
  getDistanceFn, getMetricFn} from './metric-def';

type WordInfo = {
  index: number,
  frequency: number,
};

type StemBuffer = {
  dfName: string | undefined,
  colName: string | undefined,
  dictionary: Map<string, WordInfo> | undefined,
  indices: DG.Column | undefined,
  mostCommon: Set<number> | undefined,
  metricDef: Map<string, MetricInfo> | undefined,
  aggrDistance: DISTANCE_TYPE | undefined,
};

type ColDataType = {
  type: string, 
  metric: MetricInfo, 
  use: boolean,
};

export var stemBuffer: StemBuffer = {
  dfName: undefined,
  colName: undefined,
  dictionary: undefined,
  indices: undefined,
  mostCommon: undefined,
  metricDef: undefined,
  aggrDistance: undefined,
};

const RATIO = 0.05;
const INFTY = 100000;
const TINY = 0.000001;
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
  "nor", 'nan', 'no', 'not',
  "of", "on", "once", "only", "or", "other", "ought", "our", "ours", "ourselves", "out", "over", "own", 
  "same", "she", "she'd", "she'll", "she's", "should", "so", "some", "such", 
  "than", "that", "that's", "the", "their", "theirs", "them", "themselves", "then", "there", "there's", "these", "they", "they'd", "they'll", "they're", "they've", "this", "those", "through", "to", "too", 
  "under", "until", "up", 
  "very", 
  "was", "we", "we'd", "we'll", "we're", "we've", "were", "what", "what's", "when", "when's", "where", "where's", "which", "while", "who", "who's", "whom", "why", "why's", "with", "would", 
  "you", "you'd", "you'll", "you're", "you've", "your", "yours", "yourself", "yourselves"];

const SEPARATORS = '[].,:;?!(){} \n\t#';

enum ERR_MSG {
  UNSUPPORTED_METRIC = 'Unsupported metric',
};

type FeatureSpecification = {
  source: any,
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
export function stemmColumn(source: DG.Column<string>, minCharCount: number): {dict: Map<string, WordInfo>, indices: DG.Column, mostCommon: Set<number>} {
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
      }
      else {
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

  for (const val of dict.values())
    if (val.frequency > limit)
      mostCommon.add(val.index);

  return {dict: dict, indices: indices, mostCommon: mostCommon};
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
      closestDists[i] = commonItemsCountNon(query, current) / ((queryLen < curLen) ? 1 : curLen); // */((queryLen < curLen) ? queryLen : curLen);
  }
    
  let minRelIdx = closestDists.reduce((r, v, i, a) => v >= a[r] ? r : i, -1);
  
  for (let i = count; i < size; ++i) {
    const current = indices.get(i) as Uint16Array;
    const curLen = current.length;

    if (curLen > 0) {      
      const curRelDist = commonItemsCountNon(query, current) / ((queryLen < curLen) ? 1 : curLen); // */((queryLen < curLen) ? queryLen : curLen);

      if (curRelDist > closestDists[minRelIdx]) {
        closestDists[minRelIdx] = curRelDist;
        closestInds[minRelIdx] = i;
        minRelIdx = closestDists.reduce((r, v, i, a) => v >= a[r] ? r : i, -1);
      }
    }
  }
  
  return closestInds.filter(i => i !== idx);
}

/** Return number of common elements. */
function commonItemsCountNon(arr1: Uint16Array, arr2: Uint16Array): number {
  const len1 = arr1.length;
  const len2 = arr2.length;
  let count = 0;
  let i1 = 0;
  let i2 = 0;

  while((i1 < len1) && (i2 < len2))
    if (arr1[i1] === arr2[i2]) {
      if (!stemBuffer.mostCommon?.has(arr1[i1]))
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

  const common = commonItemsCountNon(ref1, ref2);  

  if (common < 2)
  //if (common === 0)
    return INFTY;
  
  //return ((len1 < len2) ? len1 : len2) / common - 1;
  return ((len1 < len2) ? 1 : len2) / common - 1;
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
export function getMarkedStringAndCommonWordsMap(queryIdx: number, strToBeMarked: string): {marked: any[], commonWords: Map<string, number>}  {
  const size = strToBeMarked.length;

  if (size === 0)
    return {marked: [] as any[], commonWords: new Map<string, number>()
    };

  const marked = [] as any[];
  const commonWords = new Map<string, number>();

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
      const wordLow = word.toLowerCase();

      if (!STOP_WORDS.includes(wordLow)) {
        const value = stemBuffer.dictionary?.get(stemmer(wordLow));

        if (wordsOfQuery.includes(value?.index)) {
          let p = ui.inlineText([word]);
          //@ts-ignore
          $(p).css({"font-weight": "bold"});
          marked.push(p);
          const buf = word.toLocaleLowerCase();
          commonWords.set(buf, (commonWords.get(buf) ?? 0) + 1);          
        }
        else
          marked.push(word);
      }
      else
      marked.push(word);
    }

    else {
      while (SEPARATORS.includes(strToBeMarked[end])) {
        ++end;

        if (end === size)
          break;
      }

      marked.push(strToBeMarked.slice(start, end));
    }   

    start = end;
    isWord = !isWord;
  } // while

  return {marked: marked, commonWords: commonWords};
}

/** */
export function modifyMetric(df: DG.DataFrame, target: DG.Column): void {
  const colsData = new Map<string, ColDataType>();
  const colsInp = new Map<string, {metricInput: DG.InputBase, weightInput: DG.InputBase}>();
  const inputElements = new Map<string, HTMLDivElement>();
  const dlg = ui.dialog({title: 'Edit distance'});
  const initCheckedCols = [...stemBuffer.metricDef!.keys()];

  const onColumnsChanged = (columns: DG.Column[]) => {
    const names = columns.map((col) => col.name);
    
    distInput.root.hidden = (names.length < 2);
    columnsHeader.hidden = (names.length < 1);

    for (const key of colsData.keys()) {
      const val = colsData.get(key);

      if (names.includes(key)) {
        val!.use = true;
        //@ts-ignore
        inputElements.get(key)?.hidden = false;
      }
      else {
        val!.use = false;
        //@ts-ignore
        inputElements.get(key)?.hidden = true;
      }

      colsData.set(key, val!);
    }        
  };

  const colsInput = ui.columnsInput('Features', df, onColumnsChanged, {checked: initCheckedCols});  
  colsInput.setTooltip('Features used in computing similarity measure.');
  dlg.add(ui.h3('Source'));
  dlg.add(colsInput);  
  
  const distInput = ui.choiceInput('Distance', stemBuffer.aggrDistance, DISTANCE_TYPES, (dist: DISTANCE_TYPE) => { stemBuffer.aggrDistance = dist;});
  distInput.setTooltip('Type of distance between elements with the specified features.');
  dlg.add(distInput);
  distInput.root.hidden = (initCheckedCols.length < 2);

  const columnsHeader = ui.h3('Features');
  dlg.add(columnsHeader);
  columnsHeader.hidden = (initCheckedCols.length < 1);
  
  for (const col of df.columns) {
    const name = col.name;

    const colData = {
      type: col.type, 
      metric: getDefaultMetric(col),
      use: false,
    };

    if (stemBuffer.metricDef?.has(name)) {
      const val = stemBuffer.metricDef.get(name);
      colData.use = true;      
      colData.metric.type = val!.type;
      colData.metric.weight = val!.weight;
    }

    colsData.set(name, colData);

    const choices = getChoicesList(col);
    const metricInput = ui.choiceInput(`${name}:`, colData.metric.type as string, choices, (str: string) => {
      const val = colsData.get(name);
      //@ts-ignore
      val?.metric.type = str;
      colsData.set(name, val!);
    });
    metricInput.setTooltip(`Type of metric between the '${name}' feature values.`);

    const weightInput = ui.floatInput('metric with the weight', colData.metric.weight, (w: number) => {
      const val = colsData.get(name);
      //@ts-ignore
      val?.metric.weight = w;
      colsData.set(name, val!);
    });
    weightInput.setTooltip(`Weight coefficient of the '${name}' feature metric.`);
    
    const inputs = {metricInput: metricInput, weightInput: weightInput};

    colsInp.set(name, inputs);

    const uiElem = ui.divH([inputs.metricInput.root, inputs.weightInput.root]);
    inputElements.set(name, uiElem);

    uiElem.hidden = !colData.use;

    dlg.add(uiElem);
  }

  dlg  
  .onCancel(() => {})
  .onOK(() => {
    const map = new Map<string, MetricInfo>();

    for (const key of colsData.keys()) {
      const val = colsData.get(key);

      if (val?.use)
        map.set(key, val.metric);
    }

    stemBuffer.metricDef = map;
  })
  .show({x: 300, y: 300});
}

/**  */
export function setStemmingBuffer(df: DG.DataFrame, source: DG.Column): void {
  stemBuffer.dfName = df.name;
  stemBuffer.colName = source.name;

  const stemmingRes = stemmColumn(source, MIN_CHAR_COUNT);
  stemBuffer.dictionary = stemmingRes.dict;
  stemBuffer.indices = stemmingRes.indices;
  stemBuffer.mostCommon = stemmingRes.mostCommon;
  stemBuffer.metricDef = getOneHotMetricsMap(df, source);
  stemBuffer.aggrDistance = getDefaultDistnce();
}

function getTargetColDistFn(metric: MetricInfo): (query: Uint16Array, current: Uint16Array) => number {
  switch (metric.type) {
    case STR_METRIC_TYPE.STEMMING_BASED:
      return (query: Uint16Array, current: Uint16Array) => {
        const curLen = current.length;
        
        if (curLen === 0) 
          return INFTY;
        
        return metric.weight * ((query.length < curLen) ? 1 : curLen) / (commonItemsCountNon(query, current) + TINY);
      }; 
    
    case STR_METRIC_TYPE.EQUALITY:
      return (query: Uint16Array, current: Uint16Array) => {
        if (query === current)
          return metric.weight;

        return 0;
      };
    
    default:
      throw new Error(ERR_MSG.UNSUPPORTED_METRIC);
  }
}

/** 
export function getClosestAdvanced(idx: number, count: number): Uint32Array {

  console.log(stemBuffer.metricDef);
  console.log(stemBuffer.aggrDistance);  

  const indices = stemBuffer.indices!;
  const size = indices.length;
  const query = indices.get(idx) as Uint16Array;
  const queryLen = query.length;

  if ((size < count) || (queryLen === 0))
    return new Uint32Array([...Array(size).keys()]);

  const targetColMetricDef = stemBuffer.metricDef?.get(stemBuffer.colName!);  
  let targetColDistFn = (query: Uint16Array, current: Uint16Array) => 0;
  
  if (targetColMetricDef !== undefined) 
    if (targetColMetricDef.weight > 0)
      targetColDistFn = getTargetColDistFn(targetColMetricDef);

  const closestInds = new Uint32Array(count);
  const closestDists = new Float32Array(count);
    
  for (let i = 0; i < count; ++i) {
    const current = indices.get(i) as Uint16Array; 
    closestInds[i] = i;
    closestDists[i] = targetColDistFn(query, current);
  }
    
  let maxRelIdx = closestDists.reduce((r, v, i, a) => v <= a[r] ? r : i, -1);
  
  for (let i = count; i < size; ++i) {
    const current = indices.get(i) as Uint16Array;
    const curRelDist = targetColDistFn(query, current);

    if (curRelDist < closestDists[maxRelIdx]) {
      closestDists[maxRelIdx] = curRelDist;
      closestInds[maxRelIdx] = i;
      maxRelIdx = closestDists.reduce((r, v, i, a) => v <= a[r] ? r : i, -1);
    }    
  }
  
  return closestInds.filter(i => i !== idx);
}*/

/** */
export function getClosestAdvanced(df: DG.DataFrame, idx: number, count: number): Uint32Array {

  console.log(stemBuffer.metricDef);
  console.log(stemBuffer.aggrDistance);  

  const indices = stemBuffer.indices!;
  const size = indices.length;
  const queryText = indices.get(idx) as Uint16Array;
  let targetColDistFn = (a: Uint16Array, b: Uint16Array) => 0;
  let isTargetColUsed = false;
  const queryLen = queryText.length;

  if ((size < count) || (queryLen === 0) || (stemBuffer.metricDef === undefined))
    return new Uint32Array([...Array(size).keys()]);

  if (stemBuffer.metricDef.size < 1)
    return new Uint32Array([...Array(count).keys()]);

  const featuresCount = stemBuffer.metricDef.size;
  const featuresMetricVals = new Float32Array(featuresCount);
  const featureDistance = getDistanceFn(stemBuffer.aggrDistance!, featuresCount);
  const featureMap = new Map<string, FeatureSpecification>();
  const queryFeatures = new Map<string, any>();

  for (const [key, value] of stemBuffer.metricDef) 
    if (key === stemBuffer.colName!) {
      isTargetColUsed = true;
      targetColDistFn = getTargetColDistFn(value);
    }
    else {
      featureMap.set(key, {
        source: df.col(key)?.getRawData(),
        metricFn: getMetricFn(value),
      });

      queryFeatures.set(key, featureMap.get(key)?.source[idx]);
    }  

  console.log(queryFeatures);

  const closestInds = new Uint32Array(count);
  const closestDists = new Float32Array(count);
    
  for (let i = 0; i < count; ++i) {
    let j = 0;

    for (const [key, value] of featureMap) {
      featuresMetricVals[j] = value.metricFn(queryFeatures.get(key), value.source[i]);
      ++j;
    }

    if (isTargetColUsed)
      featuresMetricVals[j] = targetColDistFn(queryText, indices.get(i));

    closestInds[i] = i;
    closestDists[i] = featureDistance(featuresMetricVals);
  }
    
  let maxRelIdx = closestDists.reduce((r, v, i, a) => v <= a[r] ? r : i, -1);
  
  for (let i = count; i < size; ++i) {
    let j = 0;

    for (const [key, value] of featureMap) {
      featuresMetricVals[j] = value.metricFn(queryFeatures.get(key), value.source[i]);
      ++j;
    }

    if (isTargetColUsed)
      featuresMetricVals[j] = targetColDistFn(queryText, indices.get(i));

    const curDist = featureDistance(featuresMetricVals);

    if (curDist < closestDists[maxRelIdx]) {
      closestDists[maxRelIdx] = curDist;
      closestInds[maxRelIdx] = i;
      maxRelIdx = closestDists.reduce((r, v, i, a) => v <= a[r] ? r : i, -1);
    }    
  }
  
  return closestInds.filter(i => i !== idx);
}
