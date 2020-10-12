import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import AWS from 'aws-sdk';

export let _package = new DG.Package();

//name: supportedExt
//input: string filename
//output: bool supported
export function supportedExt(filename) {
    let extensions = `.csv, .doc, .docx, .eml, .epub, .gif, .htm, .html,
    .jpeg, .jpg, .json, .log, .mp3, .msg, .odt, .ogg, .pdf, .png, .pptx, .ps, .psv,
    .rtf, .tff, .tif, .tiff, .tsv, .txt, .wav, .xls, .xlsx`.replace(/\s{4}/g, '');
    return extensions.split(', ').some(ext => filename.endsWith(ext));
}

async function translateText(translate, params) {
    return new Promise((resolve, reject) => {
      translate.translateText(params, (err, data) => {
        if (err) {
          console.log('translateText error:', err);
          reject(err);
        }
        resolve(data.TranslatedText);
      })
    }).catch(() => {});
}

async function getCredentials() {
    let credentialsResponse = await _package.getCredentials();
    let credentials = {
      accessKeyId: credentialsResponse.parameters['accessKeyId'],
      secretAccessKey: credentialsResponse.parameters['secretAccessKey']
    };
    return credentials;
}

async function detectLanguage(textfile) {
    let langDetector = await grok.functions.eval('NLP:LanguageDetection');
    let detection = langDetector.prepare({ file: textfile });
    await detection.call();
    grok.shell.info(`Detected Language: ${detection.getParamValue('language')}`);
    return detection.getParamValue('alpha_2');
}

function testLanguagePair(sourceCode, targetCode) {
    // See the entire list of language codes at https://docs.aws.amazon.com/translate/latest/dg/what-is.html
    let supportedLanguages = ['af', 'sq', 'am', 'ar', 'az', 'bn', 'bs', 'bg', 'zh', 'zh-TW', 'hr', 'cs',
                              'da', 'fa-AF', 'nl', 'en', 'et', 'fi', 'fr', 'fr-CA', 'ka', 'de', 'el', 'ha',
                              'he', 'hi', 'hu', 'id', 'it', 'ja', 'ko', 'lv', 'ms', 'no', 'fa', 'ps', 'pl',
                              'pt', 'ro', 'ru', 'sr', 'sk', 'sl', 'so', 'es', 'es-MX', 'sw', 'sv', 'tl',
                              'ta', 'th', 'tr', 'uk', 'ur', 'vi'];
    if (!(supportedLanguages.includes(sourceCode) && supportedLanguages.includes(targetCode))) {
        return grok.shell.info('The language pair is not supported.');
    }
    if (sourceCode === targetCode) {
        return grok.shell.info('Cannot translate to the source language.');
    }
    return true;
}

async function extractText(textfile) {
    let textExtractor = await grok.functions.eval('NLP:TextExtractor');
    let extraction = textExtractor.prepare({ file: textfile });
    await extraction.call();
    let sourceText = extraction.getParamValue('text');
    // Character limit per request for real-time translation
    let maxLength = 5000;
    if (sourceText.length > maxLength) {
      grok.shell.info(`The text is too long. Translating the first ${maxLength} characters...`);
      sourceText = sourceText.substring(0, maxLength);
    }
    return sourceText;
}

//name: Translation
//tags: panel, widgets
//input: file textfile
//output: widget result
//condition: true
export async function translation(textfile) {
    // Configure AWS
    AWS.config.update({
      apiVersion: 'latest',
      credentials: await getCredentials(),
      region: 'us-east-2'
    });
    let translate = new AWS.Translate();

    // Detect the source language code and set the target language code
    let sourceCode = await detectLanguage(textfile);
    let targetCode = 'en';  // TODO: enable language suggestions in UI
    if (!testLanguagePair(sourceCode, targetCode)) return;

    // Extract text
    let sourceText = await extractText(textfile);
    if (!sourceText) return grok.shell.info('The input text is empty. Aborting...');

    grok.shell.info(`Translating the text from ${sourceCode} to ${targetCode}...`);
    let output = `Translation from ${sourceCode} to ${targetCode}:\n`;

    // Get translation
    output += await translateText(translate, {
      Text: sourceText,
      SourceLanguageCode: sourceCode,
      TargetLanguageCode: targetCode
    });

    return new DG.Widget(ui.divText(output));
}
