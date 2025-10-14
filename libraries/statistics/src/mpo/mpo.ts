/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

/// An array of [x, y] points representing the desirability line
/// [x, y] pairs are sorted by x in ascending order
export type DesirabilityLine = number[][];

/// A desirability line with its weight
export type PropertyDesirability = {
  line: DesirabilityLine;
  min?: number; /// min value of the property (optional; used for editing the line)
  max?: number; /// max value of the property (optional; used for editing the line)
  weight: number; /// 0-1
}

/// A map of desirability lines with their weights
export type DesirabilityProfile = {
  name: string;
  description: string;
  properties: { [key: string]: PropertyDesirability };
}

/// Calculates the desirability score for a given x value
/// Returns 0 if x is outside the range of the desirability line
/// Otherwise, returns the y value of the desirability line at x
export function desirabilityScore(x: number, desirabilityLine: DesirabilityLine): number {
  // If the line is empty or x is outside the range, return 0
  if (desirabilityLine.length === 0 || x < desirabilityLine[0][0] || x > desirabilityLine[desirabilityLine.length - 1][0])
    return 0;

  // Find the two points that x lies between
  for (let i = 0; i < desirabilityLine.length - 1; i++) {
    const [x1, y1] = desirabilityLine[i];
    const [x2, y2] = desirabilityLine[i + 1];

    if (x >= x1 && x <= x2) {
      // Linear interpolation between the two points
      if (x1 === x2) return y1;
      const slope = (y2 - y1) / (x2 - x1);
      return y1 + slope * (x - x1);
    }
  }

  // Should not happen if x is within bounds, but return 0 as fallback
  return 0;
}

/** Calculates the multi parameter optimization score, 0-100, 100 is the maximum */
export function mpo(dataFrame: DG.DataFrame, columns: DG.Column[]): DG.Column {
  if (columns.length === 0)
    throw new Error('No columns provided for MPO calculation.');

  const resultColumnName = dataFrame.columns.getUnusedName('MPO'); // Ensure unique name
  const resultColumn = DG.Column.float(resultColumnName, columns[0].length);

  const desirabilityTemplates = columns.map((column) => {
    const tag = column.getTag('desirabilityTemplate');
    return JSON.parse(tag) as PropertyDesirability;
  });

  resultColumn.init((i) => {
    let totalScore = 0;
    let maxScore = 0;

    for (let j = 0; j < columns.length; j++) {
      const desirability = desirabilityTemplates[j];
      const value = columns[j].get(i);
      const score = desirabilityScore(value, desirability.line);
      totalScore += desirability.weight * score;
      maxScore += desirability.weight;
    }

    return 100 * (totalScore / maxScore);
  });

  // Add the column to the table
  dataFrame.columns.add(resultColumn);
  return resultColumn;
}
