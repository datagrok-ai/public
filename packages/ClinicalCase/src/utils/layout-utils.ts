import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const STORAGE_NAME = 'clinical-case-layouts';
const MAX_STORAGE_SIZE = 5000;

/**
 * Saves layout to userStorage, splitting it into parts if it exceeds MAX_STORAGE_SIZE
 * @param studyId - Study ID
 * @param viewName - View name
 * @param layoutJson - Layout JSON string
 */
export async function saveLayoutToUserStorage(studyId: string, viewName: string, layoutJson: string): Promise<void> {
  // Delete old parts first
  await deleteLayoutFromUserStorage(studyId, viewName);

  // Calculate number of parts needed
  const totalParts = Math.ceil(layoutJson.length / MAX_STORAGE_SIZE);

  // Split and save each part in a single loop
  let partNumber = 0;
  for (let i = 0; i < layoutJson.length; i += MAX_STORAGE_SIZE) {
    const part = layoutJson.substring(i, i + MAX_STORAGE_SIZE);
    const key = `${studyId}|${viewName}|${partNumber}`;
    await grok.userSettings.add(STORAGE_NAME, key, part);
    partNumber++;
  }

  // Save metadata only if multi-part (for backward compatibility, single-part layouts don't have metadata)
  if (totalParts > 1) {
    const metadataKey = `${studyId}|${viewName}|metadata`;
    await grok.userSettings.add(STORAGE_NAME, metadataKey, JSON.stringify({parts: totalParts}));
  }
}

/**
 * Loads layout from userStorage, combining parts if necessary
 * @param studyId - Study ID
 * @param viewName - View name
 * @return Layout JSON string or null if not found
 */
export async function loadLayoutFromUserStorage(studyId: string, viewName: string): Promise<string | null> {
  try {
    // Check if metadata exists
    const metadataKey = `${studyId}|${viewName}|metadata`;
    const metadataStr = await grok.userSettings.getValue(STORAGE_NAME, metadataKey);

    if (metadataStr) {
      // Multi-part layout
      const metadata = JSON.parse(metadataStr);
      const parts: string[] = [];

      for (let partNumber = 0; partNumber < metadata.parts; partNumber++) {
        const key = `${studyId}|${viewName}|${partNumber}`;
        const part = await grok.userSettings.getValue(STORAGE_NAME, key);
        if (!part)
          return null; // Missing part

        parts.push(part);
      }

      return parts.join('');
    } else {
      // Single part layout
      const key = `${studyId}|${viewName}|0`;
      return await grok.userSettings.getValue(STORAGE_NAME, key);
    }
  } catch (e) {
    return null;
  }
}

/**
 * Deletes all parts of a saved layout from userStorage
 * @param studyId - Study ID
 * @param viewName - View name
 */
async function deleteLayoutFromUserStorage(studyId: string, viewName: string): Promise<void> {
  try {
    const metadataKey = `${studyId}|${viewName}|metadata`;
    const metadataStr = await grok.userSettings.getValue(STORAGE_NAME, metadataKey);

    if (metadataStr) {
      const metadata = JSON.parse(metadataStr);
      // Delete all parts
      for (let partNumber = 0; partNumber < metadata.parts; partNumber++) {
        const key = `${studyId}|${viewName}|${partNumber}`;
        await grok.userSettings.delete(STORAGE_NAME, key);
      }
      // Delete metadata
      await grok.userSettings.delete(STORAGE_NAME, metadataKey);
    } else {
      // Try to delete single part
      const key = `${studyId}|${viewName}|0`;
      await grok.userSettings.delete(STORAGE_NAME, key);
    }
  } catch (e) {
    // Ignore errors when deleting
  }
}

/**
 * Sets up layout save/load functionality for a table view
 * Adds a "Save layout" button to the ribbon and applies saved layout if exists
 * @param tableView - Table view to set up
 * @param studyId - Study ID
 * @param viewName - View name
 */
export async function setupTableViewLayout(tableView: DG.TableView, studyId: string, viewName: string): Promise<void> {
  const ribbons = tableView.getRibbonPanels();

  // Add save layout button
  const saveLayoutButton = ui.button('Save layout', async () => {
    try {
      const layout = tableView.saveLayout();
      const layoutJson = layout.toJson();
      await saveLayoutToUserStorage(studyId, viewName, layoutJson);
      grok.shell.info('Layout saved successfully');
    } catch (e: any) {
      grok.shell.error(`Failed to save layout: ${e?.message ?? e}`);
    }
  });
  ribbons.push([saveLayoutButton]);
  tableView.setRibbonPanels(ribbons);

  // Load and apply saved layout if exists
  try {
    const savedLayoutJson = await loadLayoutFromUserStorage(studyId, viewName);
    if (savedLayoutJson) {
      const layout = DG.ViewLayout.fromJson(savedLayoutJson);
      tableView.loadLayout(layout);
    }
  } catch (e: any) {
    // Silently fail if layout doesn't exist or can't be loaded
    console.log(`No saved layout found or failed to load: ${e?.message ?? e}`);
  }
}
