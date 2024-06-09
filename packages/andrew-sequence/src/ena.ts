import {EnaSequence} from './types';

export const parseENASequenceFasta = (data: string): EnaSequence | null => {
  // Split the input into lines
  const lines = data.split('\n').map((line) => line.trim());

  // The first line is the header
  const header = lines[0];

  // Re-join the rest of the lines as the sequence with newlines preserved
  const sequence = lines.slice(1).join('\n');

  // Regular expression to match and capture parts of the header
  const headerRegex = /^([^|]+)\|([^|]+)\|([^|]+)\s+(.+?)\s+(\S+)\s+(.+?),\s+(.+)\.$/;
  const match = header.match(headerRegex);

  if (!match) {
    console.error('Header format does not match expected pattern.');
    return null;
  }

  // Extract `code` and `description` from `codeDescription`
  const codeDescriptionParts = match[4].split(' ');
  const description = codeDescriptionParts.slice(0, -1).join(' ');
  const adjustedCode = codeDescriptionParts[codeDescriptionParts.length - 1];

  return {
    seqType: match[1],
    id: match[2],
    genBank: match[3],
    code: adjustedCode,
    description,
    name: match[6],
    extra: match[7],
    sequence: sequence,
    raw: data,
  };
};

export const searchENAIdsEmbs = (data: string): string[] => {
  // Split the content by lines
  const lines = data.split('\n');
  const enaIds: string[] = [];

  // Iterate through lines and collect ENA IDs, skipping duplicates
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    // Look for lines containing "DR   ENA;" and extract the ID
    if (line.startsWith('DR   ENA;')) {
      // Extract the ID (second part of the line after "DR   ENA;")
      const parts = line.split(';');
      if (parts.length > 2) {
        const id = parts[1].trim();
        if (i % 2 === 0) enaIds.push(id);
      }
    }
  }

  return enaIds;
};

export const parseENAMultipleSequencesFasta = (data: string): string[] => {
  const sequences: string[] = [];
  const sequenceCheck = /^[ACTGNRDactgnrd]+$/;
  const lines = data.split('\n');

  let i = -1;
  for (const line of lines) {
    if (!line.trim()) continue;
    if (sequenceCheck.test(line)) {
      sequences[i] += line;
      continue;
    }
    sequences.push('');
    ++i;
  }

  return sequences;
};
