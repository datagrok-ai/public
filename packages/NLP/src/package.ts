/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import AWS from 'aws-sdk';

import lang2code from './lang2code.json';
import code2lang from './code2lang.json';
import '../css/info-panels.css';
import {getMarkedString, setStemmingCash, getClosest, stemmColumn} from './stemming-tools/stemming-tools';
import {modifyMetric} from './stemming-tools/stemming-ui';
import {multiColReduceDimensionality} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {CLOSEST_COUNT, DELIMETER, POLAR_FREQ, TINY} from './stemming-tools/constants';
import {SentenceSimilarityViewer} from './utils/sentence-similarity-viewer';

import '../css/stemming-search.css';
import {delay} from '@datagrok-libraries/utils/src/test';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {VectorMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {DistanceAggregationMethods} from '@datagrok-libraries/ml/src/distance-matrix/types';

export const SentenceEmbeddingInvalidationTime = 8 * 60 * 60 * 1000; // 8 hours
class NLPPackage extends DG.Package {
  sentenceEmbeddingsCache: { [key: string]: Float32Array } = {};
  lastLoadTime?: number;

  constructor() {
    super();
  }

  async loadCachedEmbeddings() {
    if (this.lastLoadTime && Date.now() - this.lastLoadTime < SentenceEmbeddingInvalidationTime)
      return;
    // check the file if it exists

    if (await this.files.exists('cache/embeddingsCache.d42')) {
      const file = (await this.files.readBinaryDataFrames('cache/embeddingsCache.d42'));
      const df = file[0];
      const sentences = df.col('s')!.toList();
      const embeddings: Uint8Array[] = df.col('e')!.toList();
      // conversion to temporary Uint8Array is needed because DG returns a view on arrayBuffer that is large and undeterministic
      for (let i = 0; i < sentences.length; i++)
        this.sentenceEmbeddingsCache[sentences[i]] = new Float32Array(new Uint8Array(embeddings[i]).buffer);
    }
    this.lastLoadTime = Date.now();
  }

  async saveCachedEmbeddings() {
    const sentences = Object.keys(this.sentenceEmbeddingsCache);
    const embeddings = Object.values(this.sentenceEmbeddingsCache);
    const df = DG.DataFrame.create(sentences.length);
    df.columns.add(DG.Column.fromStrings('s', sentences));
    df.columns.add(DG.Column.fromList('byte_array', 'e', embeddings.map((e) => new Uint8Array(e.buffer, 0, e.length * 4))));
    await this.files.writeBinaryDataFrames('cache/embeddingsCache.d42', [df]);
  }
}

export const _package = new NLPPackage();

// AWS service instances
let translate: AWS.Translate;
//let comprehendMedical;

// UI components for the `Translation` panel
const sourceLangInput = ui.input.choice('', {value: 'Undetermined', items: [...Object.keys(lang2code), 'Undetermined', 'Other']});
const targetLangInput = ui.input.choice('', {value: 'English', items: [...Object.keys(lang2code), 'Choose...']});
const headerDiv = ui.divH([sourceLangInput.root, targetLangInput.root], 'nlp-header-div');
const translationArea = ui.input.textArea('', {value: ''});
translationArea.input.classList.add('nlp-translation-area');
const mainDiv = ui.divV([headerDiv, translationArea.root], 'nlp-main-div');
const mainWidget = new DG.Widget(mainDiv);
// UI components for the `Entities` panel
// let entDiv = ui.divText('{}', "nlp-entity-obj");
// let entWidget = new DG.Widget(entDiv);

let sourceLang: string; let sourceCode: string;
let sourceText: string; let cropped: boolean;

async function translateText(translate: AWS.Translate, params: { Text: string; SourceLanguageCode: any; TargetLanguageCode: any; }): Promise<{translation: string, error: number}> {
  return new Promise((resolve, reject) => {
    translate.translateText(params, (err, data) => {
      if (err) reject(err);
      resolve({translation: data.TranslatedText, error: 0});
    });
  }).catch((err) => {
    return {translation: '', error: 1};
  }) as Promise<{translation: string, error: number}>;
}

// async function detectEntities(comprehendMedical, params) {
//   return new Promise((resolve, reject) => {
//     comprehendMedical.detectEntitiesV2(params, (err, data) => {
//       if (err) reject(err);
//       // Alternatively, return the `data.Entities` array
//       resolve({entities: data, error: 0});
//     })
//   }).catch((err) => {
//     return {entities: {}, error: 1}
//   });
// }

async function getCredentials(): Promise<{accessKeyId: string, secretAccessKey: string}> {
  const credentialsResponse = await _package.getCredentials();
  if (credentialsResponse == null) {
    translationArea.value = 'Package credentials are not set.';
    // entDiv.value = 'Package credentials are not set.';
    return {accessKeyId: '', secretAccessKey: ''};
  }
  const credentials = {
    accessKeyId: credentialsResponse.parameters['accessKeyId'],
    secretAccessKey: credentialsResponse.parameters['secretAccessKey'],
  };
  return credentials;
}

async function extractText(textfile: DG.FileInfo) {
  const textExtractor = await grok.functions.eval('NLP:TextExtractor');
  const extraction = textExtractor.prepare({file: textfile, extension: textfile.extension});
  await extraction.call();
  return extraction.getParamValue('text');
}

async function detectLanguage(text: string) {
  const langDetector = await grok.functions.eval('NLP:LanguageDetector');
  const detection = langDetector.prepare({text: text});
  await detection.call();
  return [detection.getParamValue('language'),
    detection.getParamValue('alpha_2'),
    detection.getParamValue('alpha_3')];
}

function testLanguagePair(sourceCode: string, targetCode: any) {
  if (targetLangInput.value === 'Choose...') return false;
  const supportedLanguages = Object.keys(code2lang);
  if (!(supportedLanguages.includes(sourceCode))) {
    // The user unintentionally picks `Undetermined` or `Other`
    if (supportedLanguages.includes((<{[key: string]: string}>lang2code)[sourceLang])) {
      translationArea.value = `Translating from ${sourceLang}`;
      sourceLangInput.value = sourceLang;
      return true;
    }
    translationArea.value = (sourceLang === 'Undetermined') ? 'The language could not be determined.' :
      `The detected language (${sourceLang}) is not supported.`;
    return false;
  }
  if (sourceCode === targetCode) {
    targetLangInput.value = 'Choose...';
    return false;
  }
  return true;
}

async function doTranslation() {
  translationArea.value = '';
  const sourceLang = sourceLangInput.stringValue;
  const targetLang = targetLangInput.stringValue;
  const sourceCode = (<{[key: string]: string}>lang2code)[sourceLang];
  const targetCode = (<{[key: string]: string}>lang2code)[targetLang];
  if (!testLanguagePair(sourceCode, targetCode)) return;
  translationArea.value = 'Translating...';
  const output = await translateText(translate, {
    Text: sourceText,
    SourceLanguageCode: sourceCode,
    TargetLanguageCode: targetCode,
  });
  if (output.error === 1) translationArea.value = 'Error calling Amazon Translate.';
  else translationArea.value = output.translation + (cropped ? '...' : '');
}

//name: Translation
//input: file textfile
//output: widget result
//condition: isTextFile(textfile)
export async function translationPanel(textfile: DG.FileInfo) {
  sourceLangInput.onChanged.subscribe(async (_) => doTranslation());
  targetLangInput.onChanged.subscribe(async (_) => doTranslation());

  sourceText = await extractText(textfile);
  if (!sourceText) {
    sourceLangInput.value = 'Undetermined';
    translationArea.value = 'The input text is empty.';
    return mainWidget;
  }

  // Character limit per request for real-time translation
  const maxLengthBytes = 5000;
  const lengthBytes = (new TextEncoder().encode(sourceText)).length;
  if (lengthBytes > maxLengthBytes) {
    cropped = true;
    sourceText = sourceText.substring(0, Math.max(
      0, sourceText.length - (lengthBytes - maxLengthBytes)));
  }

  [sourceLang, sourceCode] = (await detectLanguage(sourceText)).slice(0, 2);
  // `Other` refers to detected languages that are not currently supported by AWS
  sourceLangInput.value = (sourceCode in code2lang) ? (<{[key: string]: string}>code2lang)[sourceCode] :
    (sourceCode === 'un') ? 'Undetermined' : 'Other';
  if ((sourceLangInput.value !== 'English') && (targetLangInput.value === 'Choose...'))
    targetLangInput.value = 'English';


  return mainWidget;
}

// name: Entities
// tags: panel, widgets
// input: file textfile
// output: widget result
// condition: isTextFile(textfile)
// export async function entitiesPanel(textfile) {

//   let text = await extractText(textfile);
//   if (!text) {
//     entDiv.innerText = 'The input text is empty.';
//     return entWidget;
//   }

//   let output = await detectEntities(comprehendMedical, {Text: text});
//   if (output.error === 1) {
//     entDiv.innerText = 'Error calling Comprehend Medical.';
//     return entWidget;
//   }

//   entDiv.innerText = JSON.stringify(output.entities, null, 2);

//   return entWidget;
// }

//name: exportFunc
//tags: init
export async function initAWS() {
  AWS.config.update({
    apiVersion: 'latest',
    credentials: await getCredentials(),
    region: 'us-east-2',
  });
  translate = new AWS.Translate();
  //comprehendMedical = new AWS.ComprehendMedical();
}

//top-menu: ML | Text Clustering...
//name: Compute Text Embeddings
//description: Compute text embeddings using UMAP
export function computeEmbds(): void {
  const table = grok.shell.t;
  const textCols = table?.columns?.bySemTypeAll('Text');
  if (!table || !textCols || textCols.length === 0) {
    grok.shell.warning('Please open table with text columns');
    return;
  }
  const colInput = ui.input.column('Text Column', {table: table, value: textCols[0], filter: (col) => col.semType === 'Text'});
  const neiInput = ui.input.int('Neighbours', {value: 25, min: 1, max: 100, nullable: false});
  const minDist = table.rowCount > 1000 ? 0.003 : 0.008;
  const minDistInput = ui.input.float('Minimum Distance', {value: minDist, nullable: false});
  ui.dialog('Cluster Text').add(ui.inputs([colInput, neiInput, minDistInput])).onOK(async () => {
    await multiColReduceDimensionality(table as any, [colInput.value as any], DimReductionMethods.UMAP,
      [VectorMetricsNames.Cosine],
      [1], [DG.Func.find({package: 'NLP', name: 'sentenceEmbeddingsPreprocessingFunction'})[0] as any],
      DistanceAggregationMethods.MANHATTAN,
      true, false, {
        dbScanEpsilon: minDistInput.value ?? 0.008, dbScanMinPts: 4, learningRate: 1, minDist: 0.3, nEpochs: 0,
        nNeighbors: neiInput.value ?? 25, preprocessingFuncArgs: [{}], randomSeed: '0', spread: 1.5, useWebGPU: true,
      }, {
        fastRowCount: 1000,
      }, DG.Func.find({package: 'Eda', name: 'dbscanPostProcessingFunction'})[0] as any, {
        epsilon: minDistInput.value ?? 0.008,
        minimumPoints: 4,
      }, VectorMetricsNames.Cosine);
  }).show();
}

//name: Stem Column
//meta.supportedSemTypes: Text
//meta.supportedDistanceFunctions: Common Items
//input: column col {semType: Text}
//input: string metric
//input: int minimumCharactersCount = 1 {defaultValue: 1; min: 0; max: 100; optional: true}
//output: object result
export function stemColumnPreprocessingFunction(col: DG.Column, metric: string, minimumCharactersCount: number) {
  const stemRes = stemmColumn(col, minimumCharactersCount);
  const entries = stemRes.indices.toList();
  const options = {mostCommon: stemRes.mostCommon};
  return {entries, options};
}

//name: Radial Coloring
//input: column col1
//input: column col2
export function radialColoring(col1: DG.Column, col2: DG.Column) {
  const df = col1.dataFrame;
  if (!df)
    return;
  const rowCount = df.rowCount;

  const markerSize = new Float32Array(rowCount);
  const markerColor = new Float32Array(rowCount);

  const xMean = col1.stats.avg;
  const xStd = col1.stats.stdev + TINY; // TINY is added to prevent division by zero
  const xRaw = col1.getRawData();
  const yMean = col2.stats.avg;
  const yStd = col2.stats.stdev + TINY;
  const yRaw = col2.getRawData();

  let xNorm: number;
  let yNorm: number;
  let radius: number;
  let angle: number;

  // Marker size & color are specified using polar coordinates
  for (let i = 0; i < rowCount; ++i) {
    // get normalized embeddings
    xNorm = (xRaw[i] - xMean) / xStd;
    yNorm = (yRaw[i] - yMean) / yStd;

    // compute polar coordinates
    radius = Math.sqrt(xNorm**2 + yNorm**2);
    angle = Math.acos(xNorm / (radius + TINY)) * (yNorm > 0 ? 1 : -1);

    // heuristics
    markerSize[i] = radius;
    markerColor[i] = Math.sin(1.0 / (TINY + Math.log(radius + TINY))) * Math.sin(POLAR_FREQ * angle);
  }

  const sizeCol = DG.Column.fromFloat32Array('embeddings size', markerSize);
  const colorCol = DG.Column.fromFloat32Array('embeddings color', markerColor);

  sizeCol.name = df.columns.getUnusedName(sizeCol.name);
  colorCol.name = df.columns.getUnusedName(colorCol.name);

  df.columns.add(sizeCol);
  df.columns.add(colorCol);
  const tv = grok.shell.tableView(df.name);
  if (!tv)
    return;
  const colNames = [col1.name, col2.name];
  for (const v of tv.viewers) {
    if (v instanceof DG.ScatterPlotViewer && colNames.includes(v.props.xColumnName) && colNames.includes(v.props.yColumnName)) {
      v.props.sizeColumnName = sizeCol.name;
      v.props.colorColumnName = colorCol.name;
      return;
    }
  }
}

//name: Distance
//input: string query {semType: Text}
//output: widget result
//condition: true
export function distance(query: string): DG.Widget {
  const df = grok.shell.t;
  const source = df.currentCol;

  setStemmingCash(df, source);

  const uiElem = ui.label('Edit');
  uiElem.classList.add('nlp-stemming-edit');
  uiElem.onclick = () => modifyMetric(df);
  ui.tooltip.bind(uiElem, 'Edit text similarity measure');
  const wgt = new DG.Widget(uiElem);

  return wgt;
}

//name: Similar
//input: string query {semType: Text}
//output: widget result
//condition: true
export function similar(query: string): DG.Widget {
  const df = grok.shell.t;
  const source = df.currentCol;
  const queryIdx = df.currentRowIdx;

  setStemmingCash(df, source);

  const closest = getClosest(df, queryIdx, CLOSEST_COUNT);

  const uiElements = [] as HTMLElement[];

  for (let i = 0; i < closest.length; ++i) {
    const uiElem = ui.inlineText(getMarkedString(closest[i], queryIdx, source.get(closest[i])));
    uiElements.push(uiElem);
    uiElements.push(ui.divText(DELIMETER));
  }

  return new DG.Widget(ui.divV(uiElements));
}

//name: Sentence Embeddings
//tags: dim-red-preprocessing-function
//meta.supportedSemTypes: Text
//meta.supportedDistanceFunctions: Vector Cosine, Euclidean, Manhattan
//input: column col {semType: Text}
//input: string metric
//output: object result
export async function sentenceEmbeddingsPreprocessingFunction(col: DG.Column, metric: string) {
  const embeddings: Float32Array[] = await getEmbeddings(col.toList());
  return {entries: embeddings, options: {}};
}


//name: getEmbeddings
//input: list<string> sentences
//output: string result
export async function getEmbeddings(sentences: string[]): Promise<Float32Array[]> {
  await _package.loadCachedEmbeddings();

  const resultingEmbeddings = new Array<Float32Array>(sentences.length).fill(null as any).map((_) => new Float32Array(384));
  const missingSentences: string[] = [];
  const missingIndices: number[] = [];
  for (let i = 0; i < sentences.length; i++) {
    if (_package.sentenceEmbeddingsCache[sentences[i]])
      resultingEmbeddings[i] = _package.sentenceEmbeddingsCache[sentences[i]];
    else {
      missingSentences.push(sentences[i]);
      missingIndices.push(i);
    }
  }

  if (missingSentences.length == 0)
    return resultingEmbeddings;

  const container = await grok.dapi.docker.dockerContainers.filter('nlp').first();

  if (!container.status.startsWith('started') && !container.status.startsWith('checking')) {
    try {
      await grok.dapi.docker.dockerContainers.run(container.id, true);
    } catch (e) {
      // not an admin
    }
    await delay(5000);
  }

  // do a warmup fetch to avoid errors
  try {
    await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/get_embeddings', {
      method: 'POST',
      body: JSON.stringify({sentences: ['some sentence']}),
      headers: {'Content-Type': 'application/json'},
    });
  } catch (e2) {
    _package.logger.error('Error warming up embeddings');
    console.error(e2);
  }

  const body = {
    sentences: missingSentences,
  };

  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/get_embeddings', {
    method: 'POST',
    body: JSON.stringify(body),
    headers: {'Content-Type': 'application/json'},
  });

  const result = await response.json();
  if (!result.embeddings)
    throw new Error('Embeddings docker returned an error');

  for (let i = 0; i < missingIndices.length; i++) {
    // update cache
    _package.sentenceEmbeddingsCache[missingSentences[i]] = new Float32Array(result.embeddings[i]);
    resultingEmbeddings[missingIndices[i]] = new Float32Array(result.embeddings[i]);
  }
  _package.saveCachedEmbeddings(); // do this lazy style in the background.

  return resultingEmbeddings;
}

export async function getEmbeddingsUsingWorker(sentences: string[]) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL('minilm-feature-worker', import.meta.url));

    worker.onmessage = (event) => {
      const {error, embedding} = event.data;
      if (embedding) {
        resolve(embedding);
        worker.terminate();
      } else if (error) {
        reject(new Error(error));
        worker.terminate();
      }
    };

    worker.onerror = (error) => {
      reject(error);
      worker.terminate();
    };

    worker.postMessage({sentences});
  });
}

//name: Sentence Similarity Search
//tags: viewer
//output: viewer result
export function sentenceSearchViewer() {
  return new SentenceSimilarityViewer();
}

//top-menu: ML | Sentence Similarity Search
//name: sentenceSearch
//output: viewer result
export function sentenceSearchTopMenu(): void {
  const view = grok.shell.tv;
  const textColumns = view.dataFrame.columns.bySemTypeAll('Text');

  const addViewer = (colName: string | null) => {
    const viewer = view.addViewer('Sentence Similarity Search', {...(colName ? {targetColumnName: colName} : {})});
    view.dockManager.dock(viewer, 'down');
  };
  if (textColumns.length === 0) {
    grok.shell.error('No text columns found.');
    return;
  }
  if (textColumns.length > 1) {
    const selector = ui.input.column('Text Column', {
      tooltipText: 'Select a text column to search for similar sentences', filter: (col) => col.semType === 'Text',
      table: view.dataFrame,
      value: textColumns[0],
    });
    ui.dialog('Select Text Column')
      .add(selector)
      .onOK(() => selector.value && addViewer(selector.value.name))
      .show();
  } else
    addViewer(textColumns[0].name);
}
