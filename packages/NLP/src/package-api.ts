import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  //Extractive Text Summarization Using TextRank Algorithm
  export async function summary(file: DG.FileInfo): Promise<string> {
    return await grok.functions.call('@datagrok/nlp:Summary', { file });
  }

  //Detect the language in which a text is written
  export async function languageDetector(text: string): Promise<string> {
    return await grok.functions.call('@datagrok/nlp:LanguageDetector', { text });
  }

  //Check how easy your text is to understand
  export async function textStatistics(file: DG.FileInfo): Promise<string> {
    return await grok.functions.call('@datagrok/nlp:TextStatistics', { file });
  }

  //Detect polar sentiments in a text
  export async function sentimentAnalysis(text: string): Promise<string> {
    return await grok.functions.call('@datagrok/nlp:SentimentAnalysis', { text });
  }

  //Sentiment classification of emotions and polarity
  export async function sentimentClassification(data: DG.DataFrame, col: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/nlp:SentimentClassification', { data, col });
  }

  //Extracts text from a file
  export async function textExtractor(file: DG.FileInfo, extension: string): Promise<string> {
    return await grok.functions.call('@datagrok/nlp:TextExtractor', { file, extension });
  }

  //Detect polar sentiments in a text with the VADER lexicon
  export async function ruleBasedSentimentAnalysis(data: DG.DataFrame, col: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/nlp:RuleBasedSentimentAnalysis', { data, col });
  }
}

export namespace funcs {
  export async function translationPanel(textfile: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:TranslationPanel', { textfile });
  }

  export async function initAWS(): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:InitAWS', {});
  }

  //Compute text embeddings using UMAP
  export async function computeEmbds(): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:ComputeEmbds', {});
  }

  export async function stemColumnPreprocessingFunction(col: DG.Column, metric: string, minimumCharactersCount: number): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:StemColumnPreprocessingFunction', { col, metric, minimumCharactersCount });
  }

  export async function radialColoring(col1: DG.Column, col2: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:RadialColoring', { col1, col2 });
  }

  export async function distance(query: string): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:Distance', { query });
  }

  export async function similar(query: string): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:Similar', { query });
  }

  export async function sentenceEmbeddingsPreprocessingFunction(col: DG.Column, metric: string): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:SentenceEmbeddingsPreprocessingFunction', { col, metric });
  }

  export async function getEmbeddings(): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:GetEmbeddings', {});
  }

  export async function sentenceSearchViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:SentenceSearchViewer', {});
  }

  export async function sentenceSearchTopMenu(): Promise<any> {
    return await grok.functions.call('@datagrok/nlp:SentenceSearchTopMenu', {});
  }
}
