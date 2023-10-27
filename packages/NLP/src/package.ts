import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import AWS from 'aws-sdk';
import lang2code from './lang2code.json';
import code2lang from './code2lang.json';
import '../css/info-panels.css';
import {stemCash, getEmbeddings, getMarkedStringAndCommonWordsMap, modifyMetric, setStemmingCash, getClosest, getEmbeddingsAdv} from './stemming-tools';

export const _package = new DG.Package();

// AWS service instances
let translate: AWS.Translate;
//let comprehendMedical;

// UI components for the `Translation` panel
const sourceLangInput = ui.choiceInput('', 'Undetermined', [...Object.keys(lang2code), 'Undetermined', 'Other']);
const targetLangInput = ui.choiceInput('', 'English', [...Object.keys(lang2code), 'Choose...']);
const headerDiv = ui.divH([sourceLangInput.root, targetLangInput.root], 'nlp-header-div');
const translationArea = ui.textInput('', '');
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
  const extraction = textExtractor.prepare({file: textfile});
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
//tags: panel, widgets
//input: file textfile
//output: widget result
//condition: isTextFile(textfile)
export async function translationPanel(textfile: DG.FileInfo) {
  sourceLangInput.onChanged(async (_: string) => doTranslation());
  targetLangInput.onChanged(async (_: string) => doTranslation());

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

//top-menu: NLP | Get dictionary...
//name: Compute Embeddings
//input: dataframe table {caption: Table; category: Data}
//input: column source {type: string; caption: Table; category: Data}
//output: dataframe dictionary
export function getDict(table: DG.DataFrame, source: DG.Column): DG.DataFrame {
  setStemmingCash(table, source);

  const dict = stemCash.dictionary;
  const size = dict?.size ?? 0;
  const words = Array<string>(size);
  const freqs = new Int32Array(size);

  let idx = 0;

  for (const item of dict!) {
    words[idx] = item[0];
    freqs[idx] = item[1].frequency;
    ++idx;
  }

  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Stemmed words', words),
    DG.Column.fromInt32Array('Frequency', freqs),
  ]);
}

//top-menu: NLP | Compute Embeddings...
//name: Compute Embeddings
//input: dataframe table {caption: Table; category: Data}
//input: column source {type: string; caption: Table; category: Data}
//input: int components = 2 {caption: Components; category: Hyperparameters} [The number of components (dimensions) to project the data to.]
//input: int epochs = 100 {caption: Epochs; category: Hyperparameters} [The number of epochs to optimize embeddings.]
//input: int neighbors = 4 {caption: Neighbors; category: Hyperparameters} [The number of nearest neighbors to construct the fuzzy manifold.]
//input: double minDist = 0.001 {caption: Minimum distance; category: Hyperparameters} [The effective minimum distance between embedded points.]
//input: double spread = 1.0 {caption: Spread; category: Hyperparameters} [The effective scale of embedded points.]
//input: bool inNewView = true {caption: New view; category: Results} [Provide results in a new view?]
//input: bool showScatter = true {caption: Scatter plot; category: Results} [Add a scatteplot with embeddings.]
export function computeEmbds(table: DG.DataFrame, source: DG.Column, components: number, epochs: number, 
  neighbors: number, minDist: number, spread: number, newView: boolean, showScatter: boolean): void 
{
  const start = new Date().getTime();
  const embds = getEmbeddingsAdv(table, source, components, epochs, neighbors, minDist, spread);
  const finish = new Date().getTime();
  console.log(`${table.name}:\n${components} components, ${epochs} epochs, ${neighbors} neibs, min_dist ${minDist}\nTime is ${finish - start}`);

  if (newView) {
    const res = DG.DataFrame.fromColumns([source, ...embds]);
    res.name = `${table.name}: ${components}D, ${epochs}E, ${neighbors}N, MD ${minDist}`;
    const view = grok.shell.addTableView(res);

    if (showScatter)
      view.addViewer(DG.VIEWER.SCATTER_PLOT);
  }
  else {
    embds.forEach(col => table.columns.add(col));
    const view = grok.shell.getTableView(table.name);
    
    if (showScatter)
      view.addViewer(DG.VIEWER.SCATTER_PLOT, {x: embds[0].name, y: embds[embds.length - 1].name});
  }
}

//top-menu: NLP | Process...
//name: Process Embeddings
//input: dataframe table
//input: column x {type: numerical}
//input: column y {type: numerical}
export function processEmds(table: DG.DataFrame, x: DG.Column, y: DG.Column) {
  const size = x.length;
  const xRaw = x.getRawData() as Float32Array;
  const yRaw = y.getRawData() as Float32Array;

  const getStats = (arr: Float32Array) => {
    let sum = 0;
    let sumOfSq = 0;
    

    for (let i = 0; i < size; ++i) {
      sum += arr[i];
      sumOfSq += arr[i] ** 2;
    }

    const mean = sum / size;

    return {mean: mean, std: Math.sqrt(sumOfSq / size - mean ** 2)};
  };

  const xStats = getStats(xRaw);
  const yStats = getStats(yRaw);

  const xMean = xStats.mean;
  const yMean = yStats.mean;

  const xStd = xStats.std;
  const yStd = yStats.std;

  const xNorm = new Float32Array(size);
  const yNorm = new Float32Array(size);
  const radius = new Float32Array(size);
  const angle = new Float32Array(size);

  const tiny = 0.000000001;

  for (let i = 0; i < size; ++i) {
    xNorm[i] = (xRaw[i] - xMean) / xStd;
    yNorm[i] = (yRaw[i] - yMean) / yStd;

    radius[i] = Math.sqrt(xNorm[i]**2 + yNorm[i]**2);
    angle[i] = Math.acos(xNorm[i] / (radius[i] + tiny)) * (yNorm[i] > 0 ? 1 : -1);
  }

  table.columns.add(DG.Column.fromFloat32Array(`${x.name}(norm)`, xNorm));
  table.columns.add(DG.Column.fromFloat32Array(`${y.name}(norm)`, yNorm));

  table.columns.add(DG.Column.fromFloat32Array(`radius`, radius));
  table.columns.add(DG.Column.fromFloat32Array(`angle`, angle));
}

//top-menu: NLP | Split...
//name: Split
//input: dataframe table
//input: column feature {type: numerical}
//input: double limit = 1.0
export function split(table: DG.DataFrame, feature: DG.Column, limit: number) {
  const featureRaw = feature.getRawData();
  const size = feature.length;
  let satisf = 0;

  featureRaw.forEach((val, idx, arr) => satisf += (val > limit ? 1 : 0));

  const nonSat = size - satisf;

  const cols = table.columns;

  const satCols = [] as DG.Column[];
  const nonSatCols = [] as DG.Column[];

  for (const c of cols) {
    satCols.push(DG.Column.fromType(c.type, c.name, satisf));
    nonSatCols.push(DG.Column.fromType(c.type, c.name, nonSat));
  }

  const satDF = DG.DataFrame.fromColumns(satCols);
  satDF.name = `${feature.name} > ${limit}`;

  const nonSatDF = DG.DataFrame.fromColumns(nonSatCols);
  nonSatDF.name = `${feature.name} <= ${limit}`;

  let satIdx = 0;
  let nonSatIdx = 0;

  for (let idx = 0; idx < size; ++idx)
    if (featureRaw[idx] > limit) {
      for (const c of cols)
        satDF.set(c.name, satIdx, c.get(idx));
        
      ++satIdx;
    }
    else {
      for (const c of cols) 
        nonSatDF.set(c.name, nonSatIdx, c.get(idx));
      
      ++nonSatIdx;
    }

  grok.shell.addTableView(satDF);
  grok.shell.addTableView(nonSatDF);
}

//name: Distance
//tags: panel, widgets
//input: string query {semType: Text}
//output: widget result
//condition: true
export function distance(query: string): DG.Widget {
  const df = grok.shell.t;
  const source = df.currentCol;

  setStemmingCash(df, source);
  
  const uiElem = ui.label('Edit');
  //@ts-ignore
  $(uiElem).css({"cursor": "pointer", "color": "#3cb173"});
  uiElem.onclick = () => { modifyMetric(df); };

  ui.tooltip.bind(uiElem, 'Edit text similarity measure.');

  const wgt = new DG.Widget(uiElem);

  return wgt;
}

//name: Similar
//tags: panel, widgets
//input: string query {semType: Text}
//output: widget result
//condition: true
export function similar(query: string): DG.Widget {
  const df = grok.shell.t;
  const source = df.currentCol;
  const queryIdx = df.currentRowIdx;

  setStemmingCash(df, source);

  const closest = getClosest(df, queryIdx, 6); 

  const uiElements = [] as HTMLElement[];
  const filterWords = new Map<string, number>();

  for (let i = 0; i < closest.length; ++i) {
    const res = getMarkedStringAndCommonWordsMap(queryIdx, source.get(closest[i]));

    for (const word of res.commonWords.keys())
      filterWords.set(word, (filterWords.get(word) ?? 0) + res.commonWords.get(word)!);

    const uiElem = ui.inlineText(res.marked);
    //@ts-ignore
    $(uiElem).css({"cursor": "pointer"});
    uiElem.onclick = () => { df.currentCell = df.cell(closest[i], source.name) };    
    uiElements.push(uiElem);    
    ui.tooltip.bind(uiElem, 'Click to navigate.');
    uiElements.push(ui.divText('________________________________'));    
  }
  
  const wgt = new DG.Widget(ui.divV(uiElements));

  return wgt;
}
