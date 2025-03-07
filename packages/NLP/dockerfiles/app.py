import logging
import numpy as np
import pandas as pd
from flask import Flask, request, jsonify
from flask_cors import CORS
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity

app = Flask(__name__)
CORS(app)

logging.basicConfig(
  level=logging.DEBUG,
  format="%(asctime)s - %(levelname)s - %(message)s",
  handlers=[
    logging.StreamHandler(),
    logging.FileHandler("app.log")
  ]
)

model = SentenceTransformer('all-MiniLM-L6-v2')

def get_sentence_embedding(text):
    return model.encode(text).tolist()

@app.route('/get_embeddings', methods=['POST'])
def get_embeddings():
    data = request.json
    sentences = data.get("sentences", [])

    if not sentences or not isinstance(sentences, list):
        return jsonify({"error": "Invalid input. Provide a list of sentences."}), 400

    sentence_embeddings = [get_sentence_embedding(s) for s in sentences]

    return jsonify({"embeddings": sentence_embeddings})

@app.route('/find_similar', methods=['POST'])
def find_similar():
    data = request.json
    embeddings = np.array(data.get("embeddings", []))
    sentences = data.get("sentences", [])
    index = data.get("index", None)

    if embeddings.shape[0] == 0 or index is None or index < 0 or index >= len(embeddings) or len(sentences) != len(embeddings):
        return jsonify({"error": "Invalid input. Provide valid embeddings, index, and sentences."}), 400

    mask = np.arange(len(sentences)) != index
    filtered_embeddings = embeddings[mask]
    filtered_sentences = [sentences[i] for i in range(len(sentences)) if i != index]

    similarities = cosine_similarity(filtered_embeddings, embeddings[index].reshape(1, -1)).flatten()

    df = pd.DataFrame({
        "sentence": filtered_sentences,
        "score": similarities,
        "idx": np.where(mask)[0]
    })

    df = df.sort_values(by="score", ascending=False)

    result = df.to_csv(index=False)

    return jsonify({"result": result})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)