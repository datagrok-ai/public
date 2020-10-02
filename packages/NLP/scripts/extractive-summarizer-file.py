#name: TextRank Summarization
#description: Extractive Text Summarization Using TextRank Algorithm
#language: python
#tags: nlp, panel
#input: file file {semType: text} [Given text]
#output: string summary {semType: text} [Summary]
#reference: https://arxiv.org/abs/1602.03606 [Variations of the Similarity Function of TextRank for Automated Summarization]
#condition: file.isFile && file.size < 1e6 && file.name.endsWith("txt")

from gensim.summarization.summarizer import summarize


with open(file) as f:
    text = " ".join(f.read().splitlines())

summary = " ".join(summarize(text, word_count=100, split=True))
