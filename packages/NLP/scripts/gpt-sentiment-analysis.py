#name: SentimentAnalysisGPT
#description: Detect polar sentiments in a text
#language: python
#tags: nlp, panel
#input: string text {semType: Text} [Statement where sentiments are expressed]
#output: string result

from openai import OpenAI
from pydantic import BaseModel, Field

class ReviewSentiment(BaseModel):
    review_text: str = Field(...)
    sentiment_score: float = Field(...)
    sentiment_label: str = Field(...)

class SentimentAnalyzer:
    def __init__(self):
        self.client = OpenAI(api_key="")
        self.model = "gpt-4o"
        self.system_prompt = (
            "You are a helpful sentiment classifier. "
            "Given a review, respond ONLY in valid JSON with these fields: "
            "review_text, sentiment_score (float from -1 to 1), and sentiment_label "
            "(one of: 'positive', 'neutral', 'negative').\n\n"
            "Be accurate when rating."
        )

    def analyze(self, review: str) -> dict:
        response = self.client.beta.chat.completions.parse(
            model=self.model,
            messages=[
                {"role": "system", "content": self.system_prompt},
                {"role": "user", "content": review},
            ],
            temperature=0.2,
            response_format=ReviewSentiment,
        )
        parsed: ReviewSentiment = response.choices[0].message.parsed
        return {
            "sentiment": parsed.sentiment_label,
            "polarity": parsed.sentiment_score
        }

analyzer = SentimentAnalyzer()
result = analyzer.analyze(text)