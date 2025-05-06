#name: SummaryGPT
#description: Extractive Text Summarization Using TextRank Algorithm
#language: python
#tags: nlp, panel
#input: file file {semType: text} [Given text]
#output: string summary {semType: text} [Summary]
#condition: file.isFile && file.size < 1e6 && supportedExt(file.name)

import textract
import cleantext
import openai
import tiktoken
from typing import List

params = {'filename': file, 'extension': 'pdf'}
text = textract.process(**params).decode().strip()
text = cleantext.clean(text, lower=False, to_ascii=False, no_line_breaks=True)

def num_tokens_from_string(string: str, model="gpt-4o") -> int:
    encoding = tiktoken.encoding_for_model(model)
    return len(encoding.encode(string))

def split_text_by_tokens(text: str, max_tokens: int = 2000, model="gpt-4o") -> List[str]:
    sentences = text.split(". ")
    chunks, current_chunk = [], ""
    for sentence in sentences:
        if num_tokens_from_string(current_chunk + sentence, model) < max_tokens:
            current_chunk += sentence + ". "
        else:
            chunks.append(current_chunk.strip())
            current_chunk = sentence + ". "
    if current_chunk:
        chunks.append(current_chunk.strip())
    return chunks

class TextSummarizer:
    def __init__(self, api_key: str, model: str = "gpt-4o", max_tokens: int = 2048):
        openai.api_key = api_key
        self.model = model
        self.max_tokens = max_tokens
        self.system_prompt = (
            "You are a helpful summarizer. "
            "Given a long text, provide a concise and accurate summary. "
            "Respond only with a summary of the text provided."
        )

    def chat_completion(self, content: str) -> str:
        messages = [
            {"role": "system", "content": self.system_prompt},
            {"role": "user", "content": f"Please summarize the following text in 80 words:\n{content}"},
        ]
        response = openai.ChatCompletion.create(
            model=self.model,
            messages=messages,
            temperature=0.3,
            max_tokens=100,
        )
        return response.choices[0].message["content"].strip()

    def summarize(self, text: str) -> str:
        chunks = split_text_by_tokens(text, max_tokens=2000, model=self.model)
        summaries = [self.chat_completion(chunk) for chunk in chunks]
        combined = " ".join(summaries)
        final_summary = self.chat_completion(combined)
        return final_summary

summarizer = TextSummarizer(api_key="")
summary = summarizer.summarize(text)
print("Final Summary:")
print(summary)