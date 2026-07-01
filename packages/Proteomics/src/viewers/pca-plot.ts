import * as DG from 'datagrok-api/dg';
import {GroupAssignment} from '../analysis/experiment-setup';
import {computePCA} from '../analysis/pca';

/** Local interface for area annotation region (may not be in installed datagrok-api types). */
interface AreaAnnotationRegionDef {
  type: string;
  area: [number, number][];
  fillColor: string;
  opacity: number;
  outlineColor: string;
  outlineWidth: number;
  x: string;
  y: string;
}

/** Group colors for PCA plot ellipses and styling. */
const GROUP_COLORS: string[] = ['#2196F3', '#FF5722'];

/** Compute 95% confidence ellipse polygon points for a set of 2D points.
 *  Returns an array of [x, y] pairs approximating the ellipse. */
function confidenceEllipse(
  xVals: number[],
  yVals: number[],
  nPoints: number = 64,
): [number, number][] | null {
  const n = xVals.length;
  if (n < 3)
    return null;

  // Compute means
  let cx = 0;
  let cy = 0;
  for (let i = 0; i < n; i++) {
    cx += xVals[i];
    cy += yVals[i];
  }
  cx /= n;
  cy /= n;

  // Compute 2x2 covariance matrix
  let sxx = 0;
  let sxy = 0;
  let syy = 0;
  for (let i = 0; i < n; i++) {
    const dx = xVals[i] - cx;
    const dy = yVals[i] - cy;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
  }
  sxx /= (n - 1);
  sxy /= (n - 1);
  syy /= (n - 1);

  // Eigendecompose 2x2 matrix analytically
  const trace = sxx + syy;
  const det = sxx * syy - sxy * sxy;
  const disc = Math.sqrt(Math.max(trace * trace / 4 - det, 0));
  const lambda1 = trace / 2 + disc;
  const lambda2 = trace / 2 - disc;

  if (lambda1 <= 0 || lambda2 <= 0)
    return null;

  // Eigenvector for lambda1
  let angle: number;
  if (Math.abs(sxy) > 1e-12)
    angle = Math.atan2(lambda1 - sxx, sxy);
  else
    angle = sxx >= syy ? 0 : Math.PI / 2;

  // Chi-squared critical value for 2 DOF at 95% confidence
  const chiSq95 = 5.991;

  // Semi-axes scaled by chi-squared
  const a = Math.sqrt(lambda1 * chiSq95);
  const b = Math.sqrt(lambda2 * chiSq95);

  // Generate ellipse polygon
  const points: [number, number][] = [];
  const cosA = Math.cos(angle);
  const sinA = Math.sin(angle);
  for (let i = 0; i < nPoints; i++) {
    const t = (2 * Math.PI * i) / nPoints;
    const ex = a * Math.cos(t);
    const ey = b * Math.sin(t);
    const x = cx + cosA * ex - sinA * ey;
    const y = cy + sinA * ex + cosA * ey;
    points.push([x, y]);
  }

  return points;
}

/** Creates a PCA plot (ScatterPlotViewer) on a sample-level DataFrame.
 *  Computes PCA from intensity columns and colors points by experimental group.
 *  Returns the viewer; the caller is responsible for adding the pcaDf to a view. */
export function createPcaPlot(
  df: DG.DataFrame,
  intensityCols: string[],
  groupAssignment: GroupAssignment,
  title?: string,
): {viewer: DG.ScatterPlotViewer; pcaDf: DG.DataFrame} {
  const {pcaDf, varianceExplained} = computePCA(df, intensityCols, groupAssignment);

  // Find PC column names (they include variance in the name).
  const pc1ColName = pcaDf.columns.toList().find((c) => c.name.startsWith('PC1'))?.name;
  const pc2ColName = pcaDf.columns.toList().find((c) => c.name.startsWith('PC2'))?.name;
  if (!pc1ColName || !pc2ColName)
    throw new Error('PCA did not produce PC1/PC2 columns — too few samples?');

  const sp = pcaDf.plot.scatter({
    x: pc1ColName,
    y: pc2ColName,
    color: 'Group',
  });

  // Label all sample points (typically <50 samples in proteomics)
  sp.props.labelColumnNames = ['SampleName'];
  sp.props.displayLabels = 'Always';

  // Attempt to add 95% confidence ellipses per group
  try {
    const pc1Col = pcaDf.col(pc1ColName)!;
    const pc2Col = pcaDf.col(pc2ColName)!;
    const groupCol = pcaDf.col('Group')!;
    const groups = [groupAssignment.group1.name, groupAssignment.group2.name];

    for (let g = 0; g < groups.length; g++) {
      const groupName = groups[g];
      const xVals: number[] = [];
      const yVals: number[] = [];
      for (let i = 0; i < pcaDf.rowCount; i++) {
        if (groupCol.get(i) === groupName) {
          xVals.push(pc1Col.get(i) as number);
          yVals.push(pc2Col.get(i) as number);
        }
      }

      const ellipsePoints = confidenceEllipse(xVals, yVals);
      if (ellipsePoints) {
        const color = GROUP_COLORS[g % GROUP_COLORS.length];
        const region: AreaAnnotationRegionDef = {
          type: 'area',
          area: ellipsePoints,
          fillColor: color,
          opacity: 0.15,
          outlineColor: color,
          outlineWidth: 1,
          x: pc1ColName,
          y: pc2ColName,
        };
        // annotationRegions API may not be available in all datagrok-api versions
        (sp.meta as any).annotationRegions?.add(region);
      }
    }
  } catch (e: any) {
    // Ellipse rendering is best-effort — the annotationRegions API may not be
    // exposed in the installed datagrok-api. Log the underlying error so this
    // doesn't silently mask unrelated bugs in the ellipse math or column reads.
    console.warn('PCA: could not add confidence ellipses:', e?.message ?? e);
  }

  if (title)
    sp.setOptions({title});

  return {viewer: sp, pcaDf};
}
