import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DesirabilityProfile, mpo} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

export async function _mpoDialog(table: DG.DataFrame): Promise<void> {
  // Get list of MPO files from appData/Chem/mpo
  const mpoFiles = await grok.dapi.files.list('System:AppData/Chem/mpo');

  const dataFrame = grok.shell.t;
  const mpoProfileEditor = new MpoProfileEditor(dataFrame);
  const defaultTemplateName = mpoFiles.length > 0 ? mpoFiles[0].fileName : null;
  // eslint-disable-next-line max-len
  const templateInput = ui.input.choice('Template', {items: mpoFiles.map((f) => f.fileName), value: defaultTemplateName!});

  // Keep track of current template data to avoid modifying original on cancel
  let currentTemplate: DesirabilityProfile | null = null;
  let currentTemplateFileName: string | null = defaultTemplateName; // Track filename

  async function loadProfile(profileFileName: string | null) {
    currentTemplateFileName = profileFileName; // Update tracked filename
    const templateFile = mpoFiles.find((f) => f.fileName === profileFileName);
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
        console.error('Failed to save template:', e);
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
      } catch (e) {
        console.error('MPO Calculation Error:', e);
        grok.shell.error(`MPO calculation failed: ${e instanceof Error ? e.message : String(e)}`);
      }
    })
    .show();

  loadProfile(templateInput.value!);
}
