import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


// TODO: For later, i guess, if the need arises, but this can be well augmented
// by importing RCSB's data-api graphql defitinitions for stuff like ligands and nonpoly entities
// and in general much smoother and precise access to data.


/**
 * GraphQL adapter for RCSB PDB API.
 * Provides efficient, selective data retrieval from RCSB PDB using GraphQL.
 */
export class RcsbGraphQLAdapter {
  private static readonly GRAPHQL_ENDPOINT = 'https://data.rcsb.org/graphql';

  private static readonly QUERIES = {
    POLYMER_SEQUENCES: `
      query GetPolymerSequences($entryId: String!) {
        entry(entry_id: $entryId) {
          rcsb_id
          polymer_entities {
            rcsb_id
            entity_poly {
              pdbx_seq_one_letter_code
              type
            }
            rcsb_polymer_entity_container_identifiers {
              entity_id
              asym_ids
            }
            rcsb_polymer_entity {
              pdbx_description
            }
          }
        }
      }
    `,

    ENTRY_INFO: `
      query GetEntryInfo($entryId: String!) {
        entry(entry_id: $entryId) {
          rcsb_id
          struct {
            title
          }
          rcsb_entry_info {
            resolution_combined
            experimental_method
          }
        }
      }
    `,

    BATCH_POLYMER_SEQUENCES: `
      query GetBatchPolymerSequences($entryIds: [String!]!) {
        entries(entry_ids: $entryIds) {
          rcsb_id
          struct {
            title
          }
          polymer_entities {
            rcsb_id
            entity_poly {
              pdbx_seq_one_letter_code
              type
            }
            rcsb_polymer_entity_container_identifiers {
              entity_id
              asym_ids
            }
            rcsb_polymer_entity {
              pdbx_description
            }
          }
        }
      }
    `
  };

  /**
   * Execute a GraphQL query against RCSB PDB API.
   */
  private static async executeQuery<T>(query: string, variables: Record<string, any>): Promise<T> {
    try {
      const response = await grok.dapi.fetchProxy(this.GRAPHQL_ENDPOINT, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          query,
          variables,
        }),
      });

      if (!response.ok)
        throw new Error(`GraphQL request failed: ${response.status} ${response.statusText}`);

      const result = await response.json();

      if (result.errors)
        throw new Error(`GraphQL errors: ${JSON.stringify(result.errors)}`);

      return result.data as T;
    } catch (error) {
      console.error('GraphQL query failed:', error);
      throw new Error(`Failed to execute GraphQL query: ${error}`);
    }
  }

  /**
   * Get polymer sequences for a PDB entry.
   */
  static async getPolymerSequences(pdbId: string): Promise<{
    [entityId: string]: {
      sequence: string;
      type: string;
      chains: string[];
    }
  }> {
    const data = await this.executeQuery<{
      entry: {
        rcsb_id: string;
        polymer_entities?: Array<{
          rcsb_id: string;
          entity_poly?: {
            pdbx_seq_one_letter_code?: string;
            type?: string;
          };
          rcsb_polymer_entity_container_identifiers?: {
            entity_id: string;
            asym_ids?: string[];
          };
          rcsb_polymer_entity?: {
            pdbx_description?: string;
          };
        }>;
      };
    }>(
      this.QUERIES.POLYMER_SEQUENCES,
      {entryId: pdbId.toUpperCase()}
    );

    if (!data.entry || !data.entry.polymer_entities)
      throw new Error(`No polymer entities found for PDB ID: ${pdbId}`);

    const result: {[entityId: string]: {sequence: string; type: string; chains: string[]}} = {};

    for (const entity of data.entry.polymer_entities) {
      if (entity.entity_poly?.pdbx_seq_one_letter_code &&
          entity.rcsb_polymer_entity_container_identifiers?.entity_id) {
        const entityId = entity.rcsb_polymer_entity_container_identifiers.entity_id;
        const chains = entity.rcsb_polymer_entity_container_identifiers.asym_ids || [];

        result[entityId] = {
          sequence: entity.entity_poly.pdbx_seq_one_letter_code.replace(/\n/g, '').trim(),
          type: entity.entity_poly.type || 'unknown',
          chains: chains
        };
      }
    }

    return result;
  }

  /**
   * Batch extract protein sequences for multiple PDB IDs.
   */
  static async batchGetProteinSequences(
    pdbIds: string[],
    batchSize: number = 10,
    onProgress?: (completed: number, total: number, currentPdbId: string) => void
  ): Promise<{[pdbId: string]: {[chainId: string]: string}}> {
    const results: {[pdbId: string]: {[chainId: string]: string}} = {};
    const uniquePdbIds = Array.from(new Set(pdbIds.map((id) => id.trim().toUpperCase())));
    let completed = 0;

    for (let i = 0; i < uniquePdbIds.length; i += batchSize) {
      const batch = uniquePdbIds.slice(i, i + batchSize);

      const batchPromises = batch.map(async (pdbId) => {
        try {
          onProgress?.(completed, uniquePdbIds.length, pdbId);
          const sequences = await this.getProteinSequences(pdbId);
          completed++;
          return {pdbId, sequences, success: true};
        } catch (error) {
          console.error(`GraphQL batch failed for ${pdbId}:`, error);
          completed++;
          return {pdbId, sequences: {}, success: false, error: error};
        }
      });

      const batchResults = await Promise.all(batchPromises);

      for (const result of batchResults)
        results[result.pdbId] = result.sequences;

      if (i + batchSize < uniquePdbIds.length)
        await new Promise((resolve) => setTimeout(resolve, 100));
    }

    return results;
  }

  /**
   * Batch extract all polymer data for multiple PDB IDs.
   */
  static async batchGetPolymerSequences(
    pdbIds: string[],
    batchSize: number = 10,
    onProgress?: (completed: number, total: number, currentPdbId: string) => void
  ): Promise<{[pdbId: string]: {[entityId: string]: {sequence: string; type: string; chains: string[]}}}> {
    const results: {[pdbId: string]: {[entityId: string]: {sequence: string; type: string; chains: string[]}}} = {};
    const uniquePdbIds = Array.from(new Set(pdbIds.map((id) => id.trim().toUpperCase())));
    let completed = 0;

    for (let i = 0; i < uniquePdbIds.length; i += batchSize) {
      const batch = uniquePdbIds.slice(i, i + batchSize);

      const batchPromises = batch.map(async (pdbId) => {
        try {
          onProgress?.(completed, uniquePdbIds.length, pdbId);
          const polymerData = await this.getPolymerSequences(pdbId);
          completed++;
          return {pdbId, polymerData, success: true};
        } catch (error) {
          console.error(`GraphQL polymer batch failed for ${pdbId}:`, error);
          completed++;
          return {pdbId, polymerData: {}, success: false, error: error};
        }
      });

      const batchResults = await Promise.all(batchPromises);

      for (const result of batchResults)
        results[result.pdbId] = result.polymerData;

      if (i + batchSize < uniquePdbIds.length)
        await new Promise((resolve) => setTimeout(resolve, 100));
    }

    return results;
  }

  /**
   * Get basic entry information for a PDB ID.
   */
  static async getEntryInfo(pdbId: string): Promise<{
    id: string;
    title?: string;
    resolution?: number;
    experimentalMethod?: string;
  }> {
    const data = await this.executeQuery<{
      entry: {
        rcsb_id: string;
        struct?: { title: string };
        rcsb_entry_info?: {
          resolution_combined?: number[];
          experimental_method?: string;
        };
      };
    }>(
      this.QUERIES.ENTRY_INFO,
      {entryId: pdbId.toUpperCase()}
    );

    if (!data.entry)
      throw new Error(`Entry not found for PDB ID: ${pdbId}`);

    return {
      id: data.entry.rcsb_id,
      title: data.entry.struct?.title,
      resolution: data.entry.rcsb_entry_info?.resolution_combined?.[0],
      experimentalMethod: data.entry.rcsb_entry_info?.experimental_method
    };
  }

  /**
   * Get protein sequences only, filtering out DNA/RNA.
   */
  static async getProteinSequences(pdbId: string): Promise<{[chainId: string]: string}> {
    const polymerData = await this.getPolymerSequences(pdbId);
    const proteinSequences: {[chainId: string]: string} = {};

    for (const [entityId, entityData] of Object.entries(polymerData)) {
      if (entityData.type?.toLowerCase().includes('protein') ||
          entityData.type === 'polypeptide(L)') {
        for (const chainId of entityData.chains)
          proteinSequences[chainId] = entityData.sequence;
      }
    }

    return proteinSequences;
  }

  /**
   * Create a DataFrame with polymer sequence data.
   */
  static async createSequenceDataFrame(pdbId: string): Promise<DG.DataFrame> {
    const polymerData = await this.getPolymerSequences(pdbId);

    const entityIds: string[] = [];
    const sequences: string[] = [];
    const types: string[] = [];
    const chainIds: string[] = [];
    const lengths: number[] = [];

    for (const [entityId, entityData] of Object.entries(polymerData)) {
      entityIds.push(entityId);
      sequences.push(entityData.sequence);
      types.push(entityData.type);
      chainIds.push(entityData.chains.join(','));
      lengths.push(entityData.sequence.length);
    }

    const df = DG.DataFrame.fromColumns([
      DG.Column.string('pdb_id', entityIds.length).init(() => pdbId.toUpperCase()),
      DG.Column.string('entity_id', entityIds.length).init((i) => entityIds[i]),
      DG.Column.string('polymer_type', types.length).init((i) => types[i]),
      DG.Column.string('chain_ids', chainIds.length).init((i) => chainIds[i]),
      DG.Column.int('sequence_length', lengths.length).init((i) => lengths[i]),
      DG.Column.string('sequence', sequences.length).init((i) => sequences[i])
    ]);

    df.getCol('sequence').semType = 'Macromolecule';
    df.getCol('sequence').setTag('alphabet', 'PT');
    df.getCol('sequence').setTag('units', 'fasta');
    df.getCol('pdb_id').semType = 'PDB_ID';

    return df;
  }
}

/**
 * Data provider function for GraphQL-based sequence retrieval.
 */
export async function getPolymerSequencesGraphQL(pdbId: string): Promise<string> {
  try {
    const df = await RcsbGraphQLAdapter.createSequenceDataFrame(pdbId);
    return df.toCsv();
  } catch (error) {
    console.error(`Failed to get sequences for ${pdbId}:`, error);
    throw error;
  }
}

/**
 * Convenience wrapper for batch extraction.
 */
export async function batchGetProteinSequencesGraphQL(
  pdbIds: string[],
  batchSize: number = 10,
  onProgress?: (completed: number, total: number, currentPdbId: string) => void
): Promise<{[pdbId: string]: {[chainId: string]: string}}> {
  return await RcsbGraphQLAdapter.batchGetProteinSequences(pdbIds, batchSize, onProgress);
}
