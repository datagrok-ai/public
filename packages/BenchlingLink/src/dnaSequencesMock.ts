// Types for DNA Sequences API
import { DnaSequence } from './dnaSequencesApi';

function randomDnaSequence(length: number): string {
  const dna = 'ATGC';
  let seq = '';
  for (let i = 0; i < length; i++)
    seq += dna[Math.floor(Math.random() * dna.length)];
  return seq;
}

export const mockDnaSequences = {
  dnaSequences: [
    ...Array.from({length: 100}, (_, i) => {
      const len = 20 + Math.floor(Math.random() * 81); // 20-100
      return {
        id: `seq_${String(i+1).padStart(3, '0')}`,
        name: `Example DNA ${i+1}`,
        bases: randomDnaSequence(len),
        isCircular: i % 2 === 0,
        length: len,
        aliases: [`DNA${i+1}`],
        annotations: [],
        apiURL: `https://benchling.com/api/v2/dna-sequences/seq_${String(i+1).padStart(3, '0')}`,
        createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
        modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
        webURL: `https://benchling.com/benchling/f/lib_55UxcIps-registry/seq_${String(i+1).padStart(3, '0')}/edit`,
      };
    })
  ],
  nextToken: undefined,
};

export const mockDnaSequence: DnaSequence = mockDnaSequences.dnaSequences[0]; 