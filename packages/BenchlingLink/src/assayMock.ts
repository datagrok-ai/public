// Types for Assay Results and Assay Runs API
import { AssayResult, AssayRun } from './assayApi';

export const mockAssayResults = {
  assayResults: [
    ...Array.from({length: 100}, (_, i) => ({
      id: `ar_${String(i+1).padStart(3, '0')}`,
      createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
      modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
      isReviewed: i % 2 === 0,
      fields: { result: 10 + i, unit: 'ng/mL' },
      schema: { id: `schema_${i+1}`, name: `ResultSchema${i+1}` },
    }))
  ],
  nextToken: undefined,
};

export const mockAssayResult: AssayResult = mockAssayResults.assayResults[0];

export const mockAssayRuns = {
  assayRuns: [
    ...Array.from({length: 100}, (_, i) => ({
      id: `run_${String(i+1).padStart(3, '0')}`,
      createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
      isReviewed: i % 2 === 0,
      fields: { runValue: 100 + i, status: i % 2 === 0 ? 'complete' : 'pending' },
      schema: { id: `schema_${i+1}`, name: `RunSchema${i+1}` },
      apiURL: `https://benchling.com/api/v2/assay-runs/run_${String(i+1).padStart(3, '0')}`,
    }))
  ],
  nextToken: undefined,
};

export const mockAssayRun: AssayRun = mockAssayRuns.assayRuns[0]; 