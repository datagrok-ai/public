#name: Summary
#description: Extractive Text Summarization Using TextRank Algorithm
#language: python
#environment: channels: [conda-forge], dependencies: [python=3.8, {pip: [gensim==3.6.0, textract, cleantext]}]
#tags: nlp, panel
#input: file file {semType: text} [Given text]
#output: string summary {semType: text} [Summary]
#reference: https://arxiv.org/abs/1602.03606 [Variations of the Similarity Function of TextRank for Automated Summarization]
#condition: file.isFile && file.size < 1e6 && supportedExt(file.name)

import textract
import cleantext
from gensim.summarization.summarizer import summarize


# Extract text
params = {'filename': file, 'extension': file[file.rfind('.', 0, -10) : -10]}
text = textract.process(**params).decode().strip()
text = cleantext.clean(text, lower=False, to_ascii=False, no_line_breaks=True)

# Return a summary or the input text if only one sentence is present
try:
    summary = " ".join(summarize(text, word_count=100, split=True))
except ValueError:
    summary = text
