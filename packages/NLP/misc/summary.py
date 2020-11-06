#name: Summarization
#language: python
#input: file file [Text file that should be summarized]
#output: string summary [Text summary]
#tags: demo, ml, nltk
#condition: file.size < 1e6 && file.name.endsWith("txt")
import pandas as pd
import numpy as np
import networkx as nx
import nltk
from nltk.corpus import stopwords
from nltk.tokenize import sent_tokenize
from nltk.cluster.util import cosine_distance

nltk.download('punkt')
nltk.download('stopwords')


def sent_similarity(x, y, stopwords=None):
    if stopwords is None:
        stopwords = []
    x = [w.lower() for w in x]
    y = [w.lower() for w in y]
    words = list(set(x + y))
    
    vec = [0] * len(words)
    vec2 = [0] * len(words)
    
    # Build vectors
    for w in x:
        if w in stopwords:
            continue
        vec[words.index(w)] += 1
    for w in y:
        if w in stopwords:
            continue
        vec2[words.index(w)] += 1
    cos_similarity = (1 - (cosine_distance(vec, vec2)))
    return cos_similarity


def make_matrix(sentences, stopwords):
    length = len(sentences)
    similarity_matrix = np.zeros((length, length))
    for z in range(length):
        for p in range(length):
            if z == p:
                continue 
            similarity_matrix[z][p] = sent_similarity(sentences[z], sentences[p])
    return similarity_matrix


def build_summary(sentences, n=5):
    stopword = stopwords.words('english')
    # Build similarity matrix
    sent_similarity_matrix = make_matrix(sentences, stopword)
    # Rank scores
    nx_graph = nx.from_numpy_array(sent_similarity_matrix)
    scores = nx.pagerank(nx_graph)
    # Sort from top to bottom
    ranked_sentence = sorted(((scores[i],s) for i, s in enumerate(sentences)), reverse=True)
    return ' '.join([ranked_sentence[i][1] for i in range(min(n, len(ranked_sentence)))])


text = open(file, 'r').read()
text = (text.encode('ascii', 'ignore')).decode("utf-8")
pd.Series(text).str.replace("[^a-zA-Z]", " ")
sentences = 3 # Number of result sentences
summary = build_summary(sent_tokenize(text), sentences)
