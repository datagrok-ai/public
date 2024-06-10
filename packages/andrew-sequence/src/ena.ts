import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface EnaSequence {
  seqType: string;
  id: string;
  genBank: string;
  sequence: string;
  code: string;
  description: string;
  name: string;
  extra: string;
  raw: string;
}

export const ENA_API = 'https://www.ebi.ac.uk/ena/browser/api/';

export function parseENASequenceFasta(data: string): EnaSequence | null {
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

export function searchENAIdsEmbs(data: string): string[] {
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

export function parseENAMultipleSequencesFasta(data: string): string[] {
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

export async function fetchSequences(query: string, limit: number, offset: number):
  Promise<{ ids: string[], sequences: string[] }> {
  const queryParams = `result=sequence&query=${query}&limit=${limit}&offset=${offset}`;
  const idsUrl = `${ENA_API}embl/textsearch?${queryParams}`;
  const seqsUrl =`${ENA_API}fasta/textsearch?${queryParams}`;

  const idsRequest: Promise<string[]> = grok.dapi.fetchProxy(idsUrl)
    .then((res: Response): Promise<string> => res.text())
    .then(searchENAIdsEmbs);

  const seqsReqest: Promise<string[]> = grok.dapi.fetchProxy(seqsUrl)
    .then((res: Response): Promise<string> => res.text())
    .then(parseENAMultipleSequencesFasta);

  const [ids, sequences] = await Promise.all([idsRequest, seqsReqest]);
  return {ids, sequences};
}

export function showEnaSequenceFormDialog(
  previewDf: DG.DataFrame, loader: (args: {
    query: string,
    limit: number,
    offset: number,
  }) => Promise<DG.DataFrame>,
) {
  const grid = DG.Viewer.grid(previewDf);
  const limitInput = ui.input.int('How many rows: ', {value: 100});
  const queryInput = ui.input.string('Query: ', {value: 'coronavirus'});
  const offsetInput = ui.input.int('Sequence offset: ', {value: 0});

  const createDf = (): Promise<DG.DataFrame> => loader({
    query: queryInput.value,
    limit: limitInput.value,
    offset: offsetInput.value,
  });

  const button = ui.button('Preview', async (): Promise<void> => {
    const df = await createDf();
    if (previewDf.rowCount > 0)
      previewDf.rows.removeAt(0, previewDf.rowCount);
    previewDf.append(df, true);
  });

  ui.dialog('Create sequences table')
    .add(ui.splitV([
      ui.splitH([
        ui.span([queryInput.root]),
        button,
      ]),
      ui.div([grid]),
      ui.div([limitInput, offsetInput]),
    ]))
    .onOK(async (): Promise<void> => {
      const df = previewDf.rowCount > 0 ? previewDf : await createDf();
      grok.shell.addTableView(df);
    })
    .show();
}
