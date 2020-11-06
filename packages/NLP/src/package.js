import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import AWS from 'aws-sdk';
import lang2code from "./lang2code.json";
import code2lang from "./code2lang.json";
import "../css/info-panels.css";

export let _package = new DG.Package();

// AWS service instances
let translate;
let comprehendMedical;

// UI components for the `Translation` panel
let sourceLangInput = ui.choiceInput('', 'Undetermined', [...Object.keys(lang2code), 'Undetermined', 'Other']);
let targetLangInput = ui.choiceInput('', 'English', [...Object.keys(lang2code), 'Choose...']);
let headerDiv = ui.div([sourceLangInput.root, ui.divText('→', "arrow"), targetLangInput.root], "header-div");
let translationArea = ui.textArea('');
translationArea.classList.add("translation-area");
let mainDiv = ui.div([headerDiv, translationArea], "main-div");
let mainWidget = new DG.Widget(mainDiv);
// UI components for the `Entities` panel
let entDiv = ui.divText('{}', "entity-obj");
let entWidget = new DG.Widget(entDiv);

let sourceLang, sourceCode;
let sourceText, cropped;

async function translateText(translate, params) {
    return new Promise((resolve, reject) => {
        translate.translateText(params, (err, data) => {
          if (err) throw err;
          resolve({ translation: data.TranslatedText, error: 0 });
        })
    }).catch((err) => { return { translation: "", error: 1 } });
}

async function detectEntities(comprehendMedical, params) {
    return new Promise((resolve, reject) => {
        comprehendMedical.detectEntitiesV2(params, (err, data) => {
            if (err) throw err;
            // Alternatively, return the `data.Entities` array
            resolve({ entities: data, error: 0 });
        })
    }).catch((err) => { return { entities: {}, error: 1 } });
}

async function getCredentials() {
    let credentialsResponse = await _package.getCredentials();
    if (credentialsResponse === null) {
        translationArea.value = 'Package credentials are not set.';
        entDiv.value = 'Package credentials are not set.';
        return {};
    }
    let credentials = {
        accessKeyId: credentialsResponse.parameters['accessKeyId'],
        secretAccessKey: credentialsResponse.parameters['secretAccessKey']
    };
    return credentials;
}

async function extractText(textfile) {
    let textExtractor = await grok.functions.eval('NLP:TextExtractor');
    let extraction = textExtractor.prepare({ file: textfile });
    await extraction.call();
    return extraction.getParamValue('text');
}

async function detectLanguage(text) {
    let langDetector = await grok.functions.eval('NLP:LanguageDetector');
    let detection = langDetector.prepare({ text: text });
    await detection.call();
    return [detection.getParamValue('language'),
            detection.getParamValue('alpha_2'),
            detection.getParamValue('alpha_3')];
}

function testLanguagePair(sourceCode, targetCode) {
    let supportedLanguages = Object.keys(code2lang);
    if (targetLangInput.value === 'Choose...') return false;
    if (!(supportedLanguages.includes(sourceCode))) {
        // The user unintentionally picks `Undetermined` or `Other`
        if (supportedLanguages.includes(lang2code[sourceLang])) {
            translationArea.value = `Translating from ${sourceLang}`;
            sourceLangInput.value = sourceLang;
            return true;
        }
        translationArea.value = `The detected language (${sourceLang}) is not supported.`;
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
    let sourceLang = sourceLangInput.stringValue;
    let targetLang = targetLangInput.stringValue;
    let sourceCode = lang2code[sourceLang];
    let targetCode = lang2code[targetLang];
    if (!testLanguagePair(sourceCode, targetCode)) return;
    translationArea.value = 'Translating...';
    let output = await translateText(translate, {
        Text: sourceText,
        SourceLanguageCode: sourceCode,
        TargetLanguageCode: targetCode
    });
    if (output.error === 1) translationArea.value = 'Error calling Amazon Translate.';
    else translationArea.value = output.translation + (cropped ? '...' : '');
}

//name: Translation
//tags: panel, widgets
//input: file textfile
//output: widget result
//condition: detectTextFile(textfile)
export async function translationPanel(textfile) {

    sourceLangInput.onChanged(async (_) => doTranslation());
    targetLangInput.onChanged(async (_) => doTranslation());
    
    sourceText = await extractText(textfile);
    if (!sourceText) {
      translationArea.value = 'The input text is empty.';
      return mainWidget;
    }

    // Character limit per request for real-time translation
    let maxLengthBytes = 5000;
    let lengthBytes = (new TextEncoder().encode(sourceText)).length;
    if (lengthBytes > maxLengthBytes) {
        cropped = true;
        sourceText = sourceText.substring(0, Math.max(
            0, sourceText.length - (lengthBytes - maxLengthBytes)));
    }

    [sourceLang, sourceCode] = (await detectLanguage(sourceText)).slice(0, 2);
    // `Other` refers to detected languages that are not currently supported by AWS
    sourceLangInput.value = sourceCode in code2lang ? code2lang[sourceCode] : 'Other';

    return mainWidget;
}

//name: Entities
//tags: panel, widgets
//input: file textfile
//output: widget result
//condition: detectTextFile(textfile)
export async function entitiesPanel(textfile) {

    let text = await extractText(textfile);
    if (!text) {
        entDiv.innerText = 'The input text is empty.';
        return entWidget;
    }

    let output = await detectEntities(comprehendMedical, { Text: text });
    if (output.error === 1) {
        entDiv.innerText = 'Error calling Comprehend Medical.';
        return entWidget;
    }

    entDiv.innerText = JSON.stringify(output.entities, null, 2);

    return entWidget;
}

//name: exportFunc
//tags: autostart
export async function initAWS() {
    AWS.config.update({
      apiVersion: 'latest',
      credentials: await getCredentials(),
      region: 'us-east-2'
    });
    translate = new AWS.Translate();
    comprehendMedical = new AWS.ComprehendMedical();
}
