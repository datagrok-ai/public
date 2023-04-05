import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import AWS from 'aws-sdk';
import lang2code from './lang2code.json';
import code2lang from './code2lang.json';
import '../css/info-panels.css';

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
