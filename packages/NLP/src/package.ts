import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import AWS from 'aws-sdk';
import lang2code from './lang2code.json';
import code2lang from './code2lang.json';
import '../css/info-panels.css';
import {getSearchResultsSplitted, getFullSearchResults, getClasses, getValidationResults,
  getStemmedColumn, getFullSearchResultsUsingStemmedData, getSearchResultsSplittedUsingStemmedData} from './stemming-tools';

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

//name: Similar
//tags: panel, widgets
//input: string query {semType: Text}
//output: widget result
//condition: true
export function similar(query: string): DG.Widget {
  const MIN_CHAR_COUNT = 1;
  const CLOSEST_COUNT = 3;

  const df = grok.shell.t;
  const baseCol = df.currentCol;

  const searchResults = getSearchResultsSplitted(query, baseCol, CLOSEST_COUNT, MIN_CHAR_COUNT);  

  const divElements = [] as HTMLDivElement[];

  for (let i = 1; i < searchResults.closestIndeces.length; ++i) {
    divElements.push(ui.divText('# ' + (searchResults.closestIndeces[i] + 1).toString()));
    divElements.push(ui.divText(searchResults.closestStrings[i]));
    divElements.push(ui.divText('------'));
  }

  const wgt = new DG.Widget(ui.divV(divElements));

  ui.tooltip.bind(wgt.root, 'Stemming-based search results.');

  return wgt;
}

//top-menu: NLP | Search...
//name: Stemming-based search
//description: Find the closest sentences with respect to stemming-based similarity metric.
//input: string query [The query string to be searched for.]
//input: dataframe table [Source dataframe.]
//input: column column [The column to be searched in.]
//input: int closest = 5 [Number of the closest items to be shown in the report.]
export function stemmingBasedSearch(query: string, table: DG.DataFrame, column: DG.Column, closest: number) { 
  const results = getFullSearchResults(query, column, closest, 1);

  results.absClosest.name = 'Search results (absolute metric)';
  results.relClosest.name = 'Search results (relative metric)';

  grok.shell.addTableView(results.absClosest);
  grok.shell.addTableView(results.relClosest);
}

//top-menu: NLP | Compare Metrics...
//name: Compare SBS-metrics
//description: Compare stemming-based similarity metrics. A query is the string specified from the target column.
//input: int index [Index of the query string taken from the target column to be searched for.]
//input: dataframe table [Source dataframe.]
//input: column column {type: string} [The column to be searched in.]
//input: column categories {type: string} [The column with validation classes.]
//input: int closest = 5 [Number of the closest items to be shown in the report.]
export function compareMetrics(index: number, table: DG.DataFrame, column: DG.Column, categories: DG.Column, closest: number) {
  const query = column.get(((index >= 0 ) && (index < column.length))? index : 0) ?? '';
  const results = getFullSearchResults(query, column, closest, 1);

  results.absClosest.name = 'Search results (absolute metric)';
  results.relClosest.name = 'Search results (relative metric)';

  results.absClosest.columns.add(getClasses(results.absClosest, categories));
  results.relClosest.columns.add(getClasses(results.relClosest, categories));

  grok.shell.addTableView(results.absClosest);
  grok.shell.addTableView(results.relClosest);
}

//top-menu: NLP | Categories Check...
//name: Stemming-based search
//description: Classes-based research of stemming-based similarity metric. A query is a random string from the column Questions.
//input: dataframe df {caption: Table; category: Database} [Source dataframe.]
//input: column questions {type: string; caption: Questions; category: Database} [The column to be searched in.]
//input: column section {type: string; caption: section; category: Database} [The column with validation classes.]
//input: int closestCount = 5 {caption: Closest count; category: Research parameters} [Number of the closest items to be shown in the report.]
//input: int launchesCount = 10 {caption: Launches; category: Research parameters} [Number of launches of the search process.]
//output: dataframe results {caption: Metrics comparison}
export function sbsMetricCheckClasses(df: DG.DataFrame, questions: DG.Column, section: DG.Column, closestCount: number, launchesCount: number) 
{
  let absClassesSum = 0;
  let relClassesSum = 0;
  let absBestSum = 0;
  let relBestSum = 0;

  for (let i = 0; i < launchesCount; ++i) {
    const idx = Math.floor(Math.random() * df.rowCount);
    const res = getValidationResults(idx, questions, section, closestCount);

    absClassesSum += res.absolute;
    relClassesSum += res.relative;
    absBestSum += res.absoluteBest ? 1 : 0;
    relBestSum += res.relativeBest ? 1 : 0;
  }

  const coef1 = 100 / launchesCount;
  const coef2 = coef1 / closestCount;

  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('criteria', [`Coincide '${section.name}', overall`, `Coincide '${section.name}', the best item`]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'absolute, %', [absClassesSum * coef2, absBestSum * coef1]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'relative, %', [relClassesSum * coef2, relBestSum * coef1]),
  ]);
}

//top-menu: NLP | Time Performance...
//name: SBS time performance
//description: Time performance of stemming-based similarity metric.
//input: dataframe table [Source dataframe.]
//input: column column {type: string} [The column to be searched in.]
//input: int index [Index of the query string taken from the target column to be searched for.]
//input: int closest = 5 [Number of the closest items to be shown in the report.]
//input: int launches = 10 [Number of the search launches.]
export function stemmingBasedSearchPerformance(table: DG.DataFrame, column: DG.Column, index: number, closest: number, launches: number) {
  const query = column.get(((index >= 0 ) && (index < column.length))? index : 0) ?? '';
  let sum = 0;

  for (let i = 0; i < launches; ++i) {
    let start = new Date().getTime();  
    const results = getFullSearchResults(query, column, closest, 1);
    let finish = new Date().getTime();

    sum += finish - start;
  }

  console.log(`Average time of stemming-based search: ${sum / launches} ms.`);

  // Q.csv, 20117 rows:
  // with removing stop words
  //   idx 587: 784.5 ms. -> 671 ms.
  //   idx 2: 799.05 ms -> 672.85 ms.
  //   idx 13367: 15603 ms. -> 901 ms.

  // without removing stop words
  //   idx 587: 761.1 ms.
  //   idx 2: 774.05 ms
  //   idx 13367: 14640 ms.
}

//top-menu: NLP | Stem column...
//name: Stem
//description: Stem column of strings and add results to the table.
//input: dataframe table
//input: column column {type: string}
export function stemColumn(table: DG.DataFrame, column: DG.Column) { table.columns.add(getStemmedColumn(column, 1)); }

//top-menu: NLP | Time Performance (using stemmed)...
//name: SBS time performance
//description: Time performance of stemming-based similarity metric.
//input: dataframe table [Source dataframe.]
//input: column column {type: string} [The column to be searched in.]
//input: int index [Index of the query string taken from the target column to be searched for.]
//input: int closest = 5 [Number of the closest items to be shown in the report.]
//input: int launches = 10 [Number of the search launches.]
export function stemmingBasedSearchPerformanceUsingStemmed(table: DG.DataFrame, column: DG.Column, index: number, closest: number, launches: number) {
  const query = column.get(((index >= 0 ) && (index < column.length))? index : 0) ?? '';
  const stemmed = table.getCol(`${column.name} (stemmed)`);
  console.log(`Looking in ${stemmed.name}`);

  let sum = 0;  

  for (let i = 0; i < launches; ++i) {
    let start = new Date().getTime();  
    const results = getFullSearchResultsUsingStemmedData(query, stemmed, closest, 1);
    let finish = new Date().getTime();

    sum += finish - start;
  }

  console.log(`Average time of stemming-based search using stemmed data: ${sum / launches} ms.`);

  // Q.csv, 20117 rows:
  // with removing stop words
  //   idx 587: 31.9 ms.
  //   idx 2: 32.6 ms.
  //   idx 13367: 219.55 ms.
}

//name: Similar (using stemmed)
//tags: panel, widgets
//input: string query {semType: Text}
//output: widget result
//condition: true
export function similarUsingStemmed(query: string): DG.Widget {
  const MIN_CHAR_COUNT = 1;
  const CLOSEST_COUNT = 3;

  const df = grok.shell.t;
  const source = df.currentCol;
  const stemmed = df.getCol(`${source.name} (stemmed)`);

  const searchResults = getSearchResultsSplittedUsingStemmedData(query, source, stemmed, CLOSEST_COUNT, MIN_CHAR_COUNT);  

  const divElements = [] as HTMLDivElement[];

  for (let i = 1; i < searchResults.closestIndeces.length; ++i) {
    divElements.push(ui.divText('# ' + (searchResults.closestIndeces[i] + 1).toString()));
    divElements.push(ui.divText(searchResults.closestStrings[i]));
    divElements.push(ui.divText('------'));
  }

  const wgt = new DG.Widget(ui.divV(divElements));

  ui.tooltip.bind(wgt.root, 'Stemming-based search results.');

  return wgt;
}