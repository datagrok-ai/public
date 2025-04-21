import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {MpoProfileEditor} from "./mpo-profile-editor";

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
export type DesirabilityProfile = {
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
      if (x1 === x2) return y1;
      const slope = (y2 - y1) / (x2 - x1);
      return y1 + slope * (x - x1);
    }
  }

  // Should not happen if x is within bounds, but return 0 as fallback
  return 0;
}

export async function _mpoDialog(table: DG.DataFrame): Promise<void> {
    // Get list of MPO files from appData/Chem/mpo
    const mpoFiles = await grok.dapi.files.list('System:AppData/Chem/mpo');

    const dataFrame = grok.shell.t;
    const mpoProfileEditor = new MpoProfileEditor(dataFrame);
    const defaultTemplateName = mpoFiles.length > 0 ? mpoFiles[0].fileName : null;
    const templateInput = ui.input.choice('Template', { items: mpoFiles.map(f => f.fileName), value: defaultTemplateName! });

    // Keep track of current template data to avoid modifying original on cancel
    let currentTemplate: DesirabilityProfile | null = null;
    let currentTemplateFileName: string | null = defaultTemplateName; // Track filename

    async function loadProfile(profileFileName: string | null) {
        currentTemplateFileName = profileFileName; // Update tracked filename
        const templateFile = mpoFiles.find(f => f.fileName === profileFileName);
        const templateContent = await templateFile!.readAsString();
        currentTemplate = JSON.parse(templateContent) as DesirabilityProfile;
        mpoProfileEditor.setProfile(currentTemplate);
    }

    // Update property info when template changes
    templateInput.onChanged.subscribe((value) => {
      loadProfile(value);
    });

    // Create dialog
    const dialog = ui.dialog('MPO Score')
        // .add(molColInput) // Molecule column is now determined automatically
        .add(templateInput)
        .add(mpoProfileEditor.root) // Add the container for properties
        .onOK(async () => {

            // --- Apply the changes ---
            // 1. Save the modified template (optional, could just use in memory)
            // Example: Save back to the file
             try {
                 const updatedTemplateString = JSON.stringify(currentTemplate, null, 2);
                 await grok.dapi.files.writeAsText(`System:AppData/Chem/mpo/${currentTemplateFileName}`, updatedTemplateString);
                 grok.shell.info(`Template '${currentTemplateFileName}' updated.`);
             } catch (e) {
                 console.error("Failed to save template:", e);
                 grok.shell.error(`Failed to save template '${currentTemplateFileName}': ${e instanceof Error ? e.message : String(e)}`);
                 // Decide if you want to proceed with calculation even if saving failed
                 // return;
             }


            // 2. Check for necessary columns and set tags
            const columnsToProcess: DG.Column[] = [];
            let missingColumns = false;
            for (const propertyName in currentTemplate!.properties) {
                const column = dataFrame.columns.byName(propertyName);
                if (!column) {
                    grok.shell.warning(`Column ${propertyName} from template not found in table. Skipping this property for calculation.`);
                    missingColumns = true;
                    continue; // Skip this property if column is missing
                }
                 // Set the *modified* desirability info as a tag
                column.setTag('desirabilityTemplate', JSON.stringify(currentTemplate!.properties[propertyName]));
                columnsToProcess.push(column);
            }

             if (columnsToProcess.length === 0) {
                 grok.shell.error('No valid columns found matching the template properties. Cannot calculate MPO score.');
                 return;
             }

             // 3. Call MPO function with the relevant columns
             try {
                const resultCol = mpo(dataFrame, columnsToProcess); // mpo now adds the column itself
                grok.shell.info(`MPO score calculated in column '${resultCol.name}'.`);
             } catch(e) {
                 console.error("MPO Calculation Error:", e);
                 grok.shell.error(`MPO calculation failed: ${e instanceof Error ? e.message : String(e)}`);
             }

        })
        .show();

    loadProfile(templateInput.value!);
}


/** Calculates the multi parameter optimization score, 0-100, 100 is the maximum */
export function mpo(dataFrame: DG.DataFrame, columns: DG.Column[]): DG.Column {
    if (columns.length === 0)
        throw new Error("No columns provided for MPO calculation.");

    const resultColumnName = dataFrame.columns.getUnusedName('MPO'); // Ensure unique name
    const resultColumn = DG.Column.float(resultColumnName, columns[0].length);

    const desirabilityTemplates = columns.map(column => {
        const tag = column.getTag('desirabilityTemplate');
        return JSON.parse(tag) as PropertyDesirability;
    });

    resultColumn.init(i => {
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