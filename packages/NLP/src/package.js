import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import AWS from 'aws-sdk';
import lang2code from "./lang2code.json";
import code2lang from "./code2lang.json";

export let _package = new DG.Package();

// UI components
let sourceLangInput = ui.choiceInput('', 'Undetermined', [...Object.keys(lang2code), 'Undetermined', 'Other']);
let targetLangInput = ui.choiceInput('', 'English', Object.keys(lang2code));
let translationArea = ui.textArea('');
let statusBar = ui.divText('');
translationArea.style = 'margin-top: 5px; padding: 5px;';
let statusBarStyle = 'color: ${1}; border: 1px solid ${2}; background: ${3}; \
                      box-shadow: 1px 1px 6px ${2}; margin-top: 5px; padding: 5px;';
let mainDiv = ui.div([
    ui.div([sourceLangInput.root, ui.divText('→'), targetLangInput.root]),
    translationArea,
    statusBar
]);
let mainWidget = new DG.Widget(mainDiv);
let isError;

function statusError(msg) {
    isError = true;
    statusBar.style = statusBarStyle.replace('${1}', '#763434')
                                    .replaceAll('${2}', '#eb6767')
                                    .replace('${3}', '#fbe0e0');
    statusBar.innerText = msg;
}

function statusReady(msg) {
    isError = false;
    statusBar.style = statusBarStyle.replace('${1}', '#286344')
                                    .replaceAll('${2}', '#3cb173')
                                    .replace('${3}', '#dcf3e7');
    statusBar.innerText = msg;
}

let sourceLang, sourceCode;
let sourceText, cropped;
let translate;

async function translateText(translate, params) {
    return new Promise((resolve, reject) => {
        translate.translateText(params, (err, data) => {
          if (err) throw err;
          resolve({ translation: data.TranslatedText, error: 0 });
        })
    }).catch((err) => { return { translation: "", error: 1 } });
}

async function getCredentials() {
    let credentialsResponse = await _package.getCredentials();
    if (credentialsResponse === null) return {};
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
    if (!(supportedLanguages.includes(sourceCode))) {
        statusError(`The detected language (${sourceLang}) is not supported.`);
        return false;
    }
    if (sourceCode === targetCode) {
        statusError('Cannot translate to the language of the original text.');
        return false;
    }
    return true;
}

async function doTranslation() {
    let sourceLang = sourceLangInput.stringValue;
    let targetLang = targetLangInput.stringValue;
    let sourceCode = lang2code[sourceLang];
    let targetCode = lang2code[targetLang];
    // Clears the text area for an unsuccessful call
    translationArea.value = '';
    if (!testLanguagePair(sourceCode, targetCode)) return;
    let output = await translateText(translate, {
        Text: sourceText,
        SourceLanguageCode: sourceCode,
        TargetLanguageCode: targetCode
    });
    if (output.error === 1) statusError('Error calling Amazon Translate.');
    else statusReady('Done!');
    translationArea.value = output.translation + (cropped ? '...' : '');
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
      statusError('The input text is empty.');
      return mainWidget;
    }

    [sourceLang, sourceCode] = (await detectLanguage(sourceText)).slice(0, 2);
    // `Other` refers to detected languages that are not currently supported by AWS
    sourceLangInput.value = sourceCode in code2lang ? code2lang[sourceCode] : 'Other';
    // Character limit per request for real-time translation
    let maxLength = 5000;
    if (sourceText.length > maxLength) {
      cropped = true;
      sourceText = sourceText.substring(0, maxLength);
    }
    doTranslation();
    return mainWidget;
}

//name: exportFunc
//tags: autostart
export async function toScriptInit() {
    // Configure AWS and create an instance of AWS Translate
    AWS.config.update({
      apiVersion: 'latest',
      credentials: await getCredentials(),
      region: 'us-east-2'
    });
    translate = new AWS.Translate();
}
