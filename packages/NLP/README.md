*NLP* is a [Datagrok package](https://datagrok.ai/help/develop/develop#packages) for natural language processing. The package provides integration with [AWS Translate](https://aws.amazon.com/translate/), a neural machine translation service, and extends Datagrok with [info panels](https://datagrok.ai/help/discover/info-panels) for text files.

*Natural Language Processing*, or *NLP* for short, is a branch of artificial intelligence that builds a bridge between computers and human languages. This field has many applications, including:

  * [language identification](https://en.wikipedia.org/wiki/Language_identification)
  * [machine translation](https://en.wikipedia.org/wiki/Machine_translation)
  * [sentiment analysis](https://en.wikipedia.org/wiki/Sentiment_analysis)
  * [text summarization](https://en.wikipedia.org/wiki/Automatic_summarization)
  * [topic modeling](https://en.wikipedia.org/wiki/Topic_model)

The package demonstrates two ways of developing [info panels](https://datagrok.ai/help/discover/info-panels) for Datagrok: with panel scripts and with JavaScript panel functions.

To write a panel script in any of the [languages supported by the platform](https://datagrok.ai/help/develop/scripting#supported-languages), you should indicate the `panel` tag and specify conditions for the panel to be shown (in the `condition` [header parameter](https://datagrok.ai/help/develop/scripting#header-parameters)):

```python
#name: Language Detection
#language: python
#input: file file {semType: text} [A text to analyze]
#output: string language {semType: lang} [Detected language]
#tags: nlp, panel
#condition: file.isFile && file.size < 1e6 && supportedExt(file.name)
```

The [scripts](https://github.com/datagrok-ai/public/tree/master/packages/NLP/scripts) folder contains more examples of such panel scripts, which are written in Python and work specifically on text files.

A different approach is used to add an info panel from a JavaScript file. The panel function should be properly annotated to return a widget. A simplified example is shown below:

```javascript
//name: Translation
//tags: panel, widgets
//input: file textfile
//output: widget result
//condition: isTextFile(textfile)
export function translationPanel(textfile) {
    return new DG.Widget(ui.divText("Lost in Translation"));
}
```

Refer to [src/package.js](https://github.com/datagrok-ai/public/blob/master/packages/NLP/src/package.js) to see the panel's complete code.

See also:

  * [Natural Language Processing](https://en.wikipedia.org/wiki/Natural_language_processing)
  * [Scripting](https://datagrok.ai/help/develop/scripting)
  * [Info Panels](https://datagrok.ai/help/discover/info-panels)
