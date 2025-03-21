import logging
import numpy as np
from flask import Flask, request, jsonify
from flask_cors import CORS
from sentence_transformers import SentenceTransformer

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

model = None

def get_sentence_embedding(text):
  if not isinstance(text, str):
    return np.zeros(384).tolist()
  return model.encode(text).tolist()

@app.route('/get_embeddings', methods=['POST'])
def get_embeddings():
  global model
  if model is None:
    logging.info("Loading model...")
    model = SentenceTransformer('all-MiniLM-L6-v2')
    logging.info("Model loaded successfully.")

  data = request.json
  sentences = data.get("sentences", [])

  if not sentences or not isinstance(sentences, list):
    return jsonify({"error": "Invalid input. Provide a list of sentences."}), 400

  sentence_embeddings = [get_sentence_embedding(s) for s in sentences]

  return jsonify({"embeddings": sentence_embeddings})

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=8000)