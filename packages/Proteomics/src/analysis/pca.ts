import * as DG from 'datagrok-api/dg';
import {GroupAssignment} from './experiment-setup';

/** Transpose a matrix (rows x cols) -> (cols x rows). */
function transpose(matrix: number[][], rows: number, cols: number): number[][] {
  const result: number[][] = [];
  for (let j = 0; j < cols; j++) {
    const row: number[] = new Array(rows);
    for (let i = 0; i < rows; i++)
      row[i] = matrix[i][j];
    result.push(row);
  }
  return result;
}

/** Multiply matrices A (m x n) and B (n x p) -> C (m x p). */
function matmul(A: number[][], B: number[][], m: number, n: number, p: number): number[][] {
  const C: number[][] = [];
  for (let i = 0; i < m; i++) {
    const row = new Array(p).fill(0);
    for (let k = 0; k < n; k++) {
      const aik = A[i][k];
      for (let j = 0; j < p; j++)
        row[j] += aik * B[k][j];
    }
    C.push(row);
  }
  return C;
}

/** Jacobi eigendecomposition for a symmetric matrix.
 *  Returns eigenvalues and eigenvectors sorted descending by eigenvalue. */
function jacobi(matrix: number[][], n: number): {eigenvalues: number[]; eigenvectors: number[][]} {
  // Deep copy the matrix
  const A: number[][] = matrix.map((row) => [...row]);

  // Initialize eigenvectors as identity matrix
  const V: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row = new Array(n).fill(0);
    row[i] = 1;
    V.push(row);
  }

  const maxIter = 100;
  const tolerance = 1e-10;
  let converged = false;

  for (let iter = 0; iter < maxIter; iter++) {
    // Find largest off-diagonal element
    let maxVal = 0;
    let p = 0;
    let q = 1;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (Math.abs(A[i][j]) > maxVal) {
          maxVal = Math.abs(A[i][j]);
          p = i;
          q = j;
        }
      }
    }

    if (maxVal < tolerance) {
      converged = true;
      break;
    }

    // Compute rotation angle
    const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
    const t = Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
    const c = 1 / Math.sqrt(t * t + 1);
    const s = t * c;

    // Apply rotation to A
    const app = A[p][p];
    const aqq = A[q][q];
    const apq = A[p][q];

    A[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
    A[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
    A[p][q] = 0;
    A[q][p] = 0;

    for (let i = 0; i < n; i++) {
      if (i !== p && i !== q) {
        const aip = A[i][p];
        const aiq = A[i][q];
        A[i][p] = c * aip - s * aiq;
        A[p][i] = A[i][p];
        A[i][q] = s * aip + c * aiq;
        A[q][i] = A[i][q];
      }
    }

    // Update eigenvectors
    for (let i = 0; i < n; i++) {
      const vip = V[i][p];
      const viq = V[i][q];
      V[i][p] = c * vip - s * viq;
      V[i][q] = s * vip + c * viq;
    }
  }

  if (!converged)
    console.warn(`PCA Jacobi did not converge in ${maxIter} iterations (n=${n}); results may be approximate`);

  // Extract eigenvalues and sort descending
  const eigenvalues = new Array(n);
  for (let i = 0; i < n; i++)
    eigenvalues[i] = A[i][i];

  const indices = Array.from({length: n}, (_, i) => i);
  indices.sort((a, b) => eigenvalues[b] - eigenvalues[a]);

  const sortedEigenvalues = indices.map((i) => eigenvalues[i]);
  const sortedEigenvectors: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row = new Array(n);
    for (let j = 0; j < n; j++)
      row[j] = V[i][indices[j]];
    sortedEigenvectors.push(row);
  }

  return {eigenvalues: sortedEigenvalues, eigenvectors: sortedEigenvectors};
}

/** Compute sample-level PCA on intensity columns.
 *  Returns a NEW sample-level DataFrame (rows = samples, not proteins)
 *  and the variance explained percentages. */
export function computePCA(
  df: DG.DataFrame,
  intensityCols: string[],
  groupAssignment: GroupAssignment,
  nComponents: number = 2,
): {pcaDf: DG.DataFrame; varianceExplained: number[]} {
  const nSamples = intensityCols.length;
  const nProteins = df.rowCount;

  // Step 1: Build transposed matrix (nSamples x nProteins)
  // Each row = one sample, each column = one protein's intensity across that sample
  const matrix: number[][] = [];

  // Compute protein means for imputation (mean across samples for each protein)
  const proteinMeans: number[] = new Array(nProteins).fill(0);
  const proteinCounts: number[] = new Array(nProteins).fill(0);

  for (const colName of intensityCols) {
    const col = df.col(colName);
    if (!col) continue;
    for (let i = 0; i < nProteins; i++) {
      if (!col.isNone(i)) {
        proteinMeans[i] += col.get(i) as number;
        proteinCounts[i]++;
      }
    }
  }
  for (let i = 0; i < nProteins; i++)
    proteinMeans[i] = proteinCounts[i] > 0 ? proteinMeans[i] / proteinCounts[i] : 0;

  // Build sample rows with imputation
  for (const colName of intensityCols) {
    const col = df.col(colName);
    const row: number[] = new Array(nProteins);
    for (let i = 0; i < nProteins; i++) {
      if (!col || col.isNone(i))
        row[i] = proteinMeans[i];
      else
        row[i] = col.get(i) as number;
    }
    matrix.push(row);
  }

  // Step 2: Center each column (protein) by subtracting its mean across samples
  const colMeans: number[] = new Array(nProteins).fill(0);
  for (let j = 0; j < nProteins; j++) {
    for (let i = 0; i < nSamples; i++)
      colMeans[j] += matrix[i][j];
    colMeans[j] /= nSamples;
  }
  for (let i = 0; i < nSamples; i++) {
    for (let j = 0; j < nProteins; j++)
      matrix[i][j] -= colMeans[j];
  }

  // Step 3: Compute covariance matrix S (nSamples x nSamples)
  // S = X * X^T / (nProteins - 1) -- small matrix since nSamples << nProteins
  const Xt = transpose(matrix, nSamples, nProteins);
  const S = matmul(matrix, Xt, nSamples, nProteins, nSamples);
  const divisor = nProteins > 1 ? nProteins - 1 : 1;
  for (let i = 0; i < nSamples; i++) {
    for (let j = 0; j < nSamples; j++)
      S[i][j] /= divisor;
  }

  // Step 4: Eigendecompose S
  const {eigenvalues, eigenvectors} = jacobi(S, nSamples);

  // Step 5: PC scores are the eigenvectors scaled by sqrt(eigenvalue)
  const scores: number[][] = [];
  for (let i = 0; i < nSamples; i++) {
    const row: number[] = [];
    for (let k = 0; k < Math.min(nComponents, nSamples); k++) {
      const scale = eigenvalues[k] > 0 ? Math.sqrt(eigenvalues[k]) : 0;
      row.push(eigenvectors[i][k] * scale);
    }
    scores.push(row);
  }

  // Step 6: Compute variance explained. Clamp eigenvalues to >= 0 for both
  // totalVariance and per-PC pct so numerical noise can't produce negative percentages.
  const totalVariance = eigenvalues.reduce((sum, v) => sum + Math.max(v, 0), 0);
  const varianceExplained: number[] = [];
  for (let k = 0; k < Math.min(nComponents, nSamples); k++) {
    const pct = totalVariance > 0 ? (Math.max(eigenvalues[k], 0) / totalVariance) * 100 : 0;
    varianceExplained.push(Math.round(pct * 10) / 10);
  }

  // Step 7: Create sample-level DataFrame. PCA on <2 samples is meaningless;
  // produce stable column labels rather than "PC1 (NaN%)".
  const pc1Pct = varianceExplained.length > 0 && isFinite(varianceExplained[0]) ? varianceExplained[0] : 0;
  const pc1Name = nSamples >= 2 ? `PC1 (${pc1Pct}%)` : 'PC1';
  const pc2Pct = varianceExplained.length >= 2 && isFinite(varianceExplained[1]) ? varianceExplained[1] : 0;
  const pc2Name = nComponents >= 2 && varianceExplained.length >= 2 && nSamples >= 2
    ? `PC2 (${pc2Pct}%)`
    : 'PC2';

  const sampleNames = DG.Column.fromStrings('SampleName', intensityCols);
  const pc1Col = DG.Column.fromFloat32Array(pc1Name,
    new Float32Array(scores.map((s) => s[0])));
  const columns: DG.Column[] = [sampleNames, pc1Col];

  if (nComponents >= 2) {
    const pc2Col = DG.Column.fromFloat32Array(pc2Name,
      new Float32Array(scores.map((s) => s.length > 1 ? s[1] : 0)));
    columns.push(pc2Col);
  }

  // Map each sample to its group
  const groupNames: string[] = intensityCols.map((colName) => {
    if (groupAssignment.group1.columns.includes(colName))
      return groupAssignment.group1.name;
    if (groupAssignment.group2.columns.includes(colName))
      return groupAssignment.group2.name;
    return 'Unknown';
  });
  const groupCol = DG.Column.fromStrings('Group', groupNames);
  columns.push(groupCol);

  const pcaDf = DG.DataFrame.fromColumns(columns);
  pcaDf.name = 'PCA';

  return {pcaDf, varianceExplained};
}
