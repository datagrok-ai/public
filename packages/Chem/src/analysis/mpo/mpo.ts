import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

/// An array of [x, y] points representing the desirability line
/// [x, y] pairs are sorted by x in ascending order
type DesirabilityLine = number[][];

/// A desirability line with its weight
export type PropertyDesirability = {
  line: DesirabilityLine;
  min?: number;      /// min value of the property (optional; used for editing the line)
  max?: number;      /// max value of the property (optional; used for editing the line)
  weight: number;   /// 0-1
}

/// A map of desirability lines with their weights
export type DesirabilityTemplate = {
  name: string;
  description: string;
  properties: { [key: string]: PropertyDesirability };
}

/// Calculates the desirability score for a given x value
/// Returns 0 if x is outside the range of the desirability line
/// Otherwise, returns the y value of the desirability line at x
function desirabilityScore(x: number, desirabilityLine: DesirabilityLine): number {
  // If the line is empty or x is outside the range, return 0
  if (desirabilityLine.length === 0 || x < desirabilityLine[0][0] || x > desirabilityLine[desirabilityLine.length - 1][0]) 
      return 0;

  // Find the two points that x lies between
  for (let i = 0; i < desirabilityLine.length - 1; i++) {
    const [x1, y1] = desirabilityLine[i];
    const [x2, y2] = desirabilityLine[i + 1];

    if (x >= x1 && x <= x2) {
      // Linear interpolation between the two points
      const slope = (y2 - y1) / (x2 - x1);
      return y1 + slope * (x - x1);
    }
  }

  return 0;
}


export function mpo(dataFrame: DG.DataFrame, columns: DG.Column[]): DG.Column {
  const resultColumn = DG.Column.float('mpo', columns[0].length);
  const desirabilityTemplates = columns.map(column => {
    return JSON.parse(column.getTag('desirabilityTemplate')) as PropertyDesirability;
  });

  resultColumn.init(i => {
    let sum = 0;
    let totalWeight = 0;

    for (let j = 0; j < columns.length; j++) {
      const value = columns[j].get(i);
      if (value !== null) {
        const desirability = desirabilityTemplates[j];
        const score = desirabilityScore(value, desirability.line);
        sum += score * desirability.weight;
        totalWeight += desirability.weight;
      }
    }

        // Normalize by total weight to get a score between 0 and 1
    return totalWeight > 0 ? sum / totalWeight : 0;
  });

  dataFrame.columns.add(resultColumn);
  return resultColumn;
}


export async function _mpoDialog(table: DG.DataFrame): Promise<void> {
  // Get list of MPO files from appData/Chem/mpo
  const mpoFiles = await grok.dapi.files.list('System:AppData/Chem/mpo');

  const dataFrame = grok.shell.t;
  const molCol = dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (!molCol)
    throw new Error('No molecule column found');

  const molColInput = ui.input.column('Column', { table: dataFrame, value: molCol })
  const templateInput = ui.input.choice('MPO Template', { items: mpoFiles.map(f => f.fileName), value: mpoFiles[0]?.fileName });
  const propertyInfos = ui.divV([]);

  // Create dialog
  const dialog = ui.dialog('MPO Score')
    .add(molColInput)
    .add(templateInput)
    .add(propertyInfos)
    .onOK(async () => {
      const templateContent = await grok.dapi.files.readAsText('System:AppData/Chem/mpo/' + templateInput.value!);
      const template = JSON.parse(templateContent) as DesirabilityTemplate;

      for (const property in template.properties) {
        const column = dataFrame.columns.byName(property);
        if (!column)
          throw new Error(`Column ${property} not found`);
        column.setTag('desirabilityTemplate', JSON.stringify(template.properties[property]));
      }

      mpo(dataFrame, Object.keys(template.properties).map(property => dataFrame.columns.byName(property)));
    })
    .show();

      // Function to update property information
  async function updatePropertyInfo(templateFileName: string) {
    const templateFile = mpoFiles.find(f => f.fileName === templateFileName);
    if (!templateFile) return;

    const templateContent = await templateFile.readAsString();
    const template = JSON.parse(templateContent) as DesirabilityTemplate;
    console.log(template);

    const grid = ui.table(Object.keys(template.properties), (propertyName, index) => {
      const prop = template.properties[propertyName];
      
      // Create canvas for the desirability line
      const canvas = ui.canvas(200, 50);
      const ctx = canvas.getContext('2d');
      if (!ctx) return [propertyName, prop.weight, canvas];
      
      // Set up canvas
      ctx.strokeStyle = '#000';
      ctx.lineWidth = 2;
      
      // Calculate scaling factors
      const minX = prop.min ?? Math.min(...prop.line.map(p => p[0]));
      const maxX = prop.max ?? Math.max(...prop.line.map(p => p[0]));
      const scaleX = 200 / (maxX - minX);
      const scaleY = 50; // Since y ranges from 0 to 1
      
      // Draw the line
      ctx.beginPath();
      prop.line.forEach((point, i) => {
        const x = (point[0] - minX) * scaleX;
        const y = 50 - (point[1] * scaleY); // Flip y-axis since canvas y increases downward
        if (i === 0) {
          ctx.moveTo(x, y);
        } else {
          ctx.lineTo(x, y);
        }
      });
      ctx.stroke();
      
      return [propertyName, prop.weight, canvas];
    })

    ui.empty(propertyInfos);
    propertyInfos.append(grid);
  }

  // Update property info when template changes
  templateInput.onChanged.subscribe((value) => {
    if (value)
      updatePropertyInfo(value);
  });
    
}
