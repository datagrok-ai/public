import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RcsbGraphQLAdapter} from './rcsb-gql-adapter';

export async function extractSequenceColumns(pdbIdColumn: DG.Column): Promise<void> {
  if (!pdbIdColumn) {
    grok.shell.error('Function was called without a valid PDB column.');
    return;
  }
  const dataFrame = pdbIdColumn.dataFrame;


  const pi = DG.TaskBarProgressIndicator.create('Extracting protein sequences via GraphQL...');

  try {
    const uniquePdbIds = pdbIdColumn.categories
      .filter((id) => id && typeof id === 'string' && id.trim().length >= 4)
      .map((id) => id.trim().toUpperCase());

    if (uniquePdbIds.length === 0) {
      grok.shell.warning('No valid PDB IDs found in the specified column.');
      return;
    }

    grok.shell.info(`Processing ${uniquePdbIds.length} unique PDB ID(s) with GraphQL...`);

    const allSequenceData = await RcsbGraphQLAdapter.batchGetProteinSequences(
      uniquePdbIds,
      20,
      (completed, total, currentPdbId) => {
        const progress = (completed / total) * 95;
        pi.update(progress, `Extracting: ${currentPdbId} (${completed}/${total})`);
      }
    );

    let maxChainCount = 0;
    for (const pdbId in allSequenceData) {
      const sequences = allSequenceData[pdbId];
      if (sequences)
        maxChainCount = Math.max(maxChainCount, Object.keys(sequences).length);
    }

    if (maxChainCount === 0) {
      grok.shell.warning('No protein sequences could be extracted from any PDB structure.');
      return;
    }

    const sequenceColumns: DG.Column[] = [];
    for (let i = 0; i < maxChainCount; i++) {
      const baseName = `Chain ${i + 1}`;
      const colName = dataFrame.columns.getUnusedName(baseName);
      sequenceColumns.push(DG.Column.string(colName, dataFrame.rowCount));
    }

    for (let i = 0; i < dataFrame.rowCount; i++) {
      const pdbId = pdbIdColumn.get(i);
      if (pdbId && allSequenceData[pdbId]) {
        const sequences = allSequenceData[pdbId];
        const chainIds = Object.keys(sequences).sort();
        for (let j = 0; j < chainIds.length && j < sequenceColumns.length; j++)
          sequenceColumns[j].set(i, sequences[chainIds[j]]);
      }
    }

    for (const column of sequenceColumns) {
      column.semType = 'Macromolecule';
      column.setTag('alphabet', 'PT');
      column.setTag('units', 'fasta');
      column.setTag(DG.TAGS.CELL_RENDERER, 'sequence');
      dataFrame.columns.add(column);
    }
  } catch (err) {
    grok.shell.error('Failed to fetch protein sequences.');
    console.error(err);
  } finally {
    pi.close();
  }
}
