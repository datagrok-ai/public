NLP is a [Datagrok package](https://datagrok.ai/help/develop/develop#packages) for natural language processing.

**Natural Language Processing**, or **NLP** for short, is a branch of artificial intelligence that builds a bridge between computers and human languages. This field has many applications, including:

- [language identification](https://github.com/datagrok-ai/public/blob/master/packages/NLP/scripts/language-detection.py)
- [machine translation](https://github.com/datagrok-ai/public/blob/master/packages/NLP/src/package.js)
- [sentiment analysis](https://github.com/datagrok-ai/public/blob/master/packages/NLP/scripts/vader-sentiment-analysis.py)
- [text summarization](https://github.com/datagrok-ai/public/blob/master/packages/NLP/scripts/extractive-summarizer.py)
- topic modeling

The package demonstrates two main ways of developing [info panels](https://datagrok.ai/help/discover/info-panels) for Datagrok: with panel scripts and with JavaScript panel functions.

To write a panel script in any of the [languages supported by the platform](https://datagrok.ai/help/develop/scripting#supported-languages), you would need to indicate the `panel` tag and specify conditions for the panel to be shown in the `condition` header parameter. For example, the [scripts](https://github.com/datagrok-ai/public/tree/master/packages/NLP/scripts) folder contains panel scripts written in Python that work on textfiles specifically.

A different approach is used to add an info panel from a JavaScript file. The panel function should be properly annotated to return a widget. A simplified example is shown below:

```javascript
//name: Translation
//tags: panel, widgets
//input: file textfile
//output: widget result
//condition: detectTextFile(textfile)
export function translationPanel(textfile) {
    return new DG.Widget(ui.divText("Lost in Translation"));
}
```

Refer to [src/package.js](https://github.com/datagrok-ai/public/blob/master/packages/NLP/src/package.js) to see the panel's complete code.
