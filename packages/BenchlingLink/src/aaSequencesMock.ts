// Types for AA Sequences API
import { AaSequence, AaSequencesPaginatedList } from './aaSequencesApi';

function randomAminoAcidSequence(length: number): string {
  const aa = 'ACDEFGHIKLMNPQRSTVWY';
  let seq = '';
  for (let i = 0; i < length; i++)
    seq += aa[Math.floor(Math.random() * aa.length)];
  return seq;
}

export const mockAaSequences: AaSequencesPaginatedList = {
  aaSequences: [
    ...Array.from({length: 100}, (_, i) => {
      const len = 20 + Math.floor(Math.random() * 81); // 20-100
      return {
        id: `prtn_${String(i+1).padStart(3, '0')}`,
        name: `Example Protein ${i+1}`,
        aminoAcids: randomAminoAcidSequence(len),
        length: len,
        aliases: [`Prot${i+1}`],
        annotations: [],
        apiURL: `https://benchling.com/api/v2/aa-sequences/prtn_${String(i+1).padStart(3, '0')}`,
        createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
        modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
        creator: {
          handle: `user${i+1}`,
          id: `ent_${i+1}`,
          name: `User ${i+1}`,
        },
        webURL: `https://benchling.com/benchling/f/lib_55UxcIps-registry/prtn_${String(i+1).padStart(3, '0')}/edit`,
      };
    })
  ],
  nextToken: undefined,
};

export const mockAaSequence: AaSequence = mockAaSequences.aaSequences[0]; 