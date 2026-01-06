/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import OpenAI from 'openai';
import {OpenAIClient} from '../openAI-client';

// no longer in use, always use core embeddings instead. keeping code for reference

export async function embedConnectionQueries(connectionId: string) {
  const queries = (await grok.dapi.queries.filter(`connection.id = "${connectionId}"`).list())
    .filter((q) => q.query?.trim());
  const openai = OpenAIClient.getInstance().openai;
  const output: {query: string, embedding: number[]}[] = [];
  const pg = DG.TaskBarProgressIndicator.create('Generating embeddings for queries...');
  let count = 0;
  for (const query of queries) {
    const emb = await getVectorEmbedding(openai, query.query!);
    output.push({query: query.query!, embedding: emb});
    count++;
    pg.update((count / queries.length) * 100, 'Generating embeddings for queries...');
  }
  pg.close();
  return output;
}

export async function getVectorEmbedding(openai: OpenAI, text: string): Promise<number[]> {
  const response = await openai.embeddings.create({
    model: 'text-embedding-3-small',
    input: text,
  });
  return response.data[0].embedding;
}

export function cosineDistance(vecA: number[], vecB: number[]): number {
  let dotProduct = 0;
  let normA = 0;
  let normB = 0;
  for (let i = 0; i < vecA.length; i++) {
    dotProduct += vecA[i] * vecB[i];
    normA += vecA[i] * vecA[i];
    normB += vecB[i] * vecB[i];
  }
  return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
}

export function getTopKSimilarQueries(textEmbedding: number[],
  candidates: {query: string, embedding: number[]}[], topK: number, threshold = 0.3): {query: string, score: number}[] {
  const scores: {query: string, score: number}[] = candidates.map((candidate) => {
    return {
      query: candidate.query,
      score: cosineDistance(textEmbedding, candidate.embedding),
    };
  });
  const out = scores.filter((s) => s.score > threshold);
  out.sort((a, b) => b.score - a.score);
  return out.slice(0, topK);
}
