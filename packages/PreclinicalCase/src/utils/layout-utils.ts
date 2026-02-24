import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const STORAGE_NAME = 'preclinical-case-layouts';
const MAX_STORAGE_SIZE = 5000;

export async function saveLayoutToUserStorage(studyId: string, viewName: string, layoutJson: string): Promise<void> {
  await deleteLayoutFromUserStorage(studyId, viewName);

  const totalParts = Math.ceil(layoutJson.length / MAX_STORAGE_SIZE);

  let partNumber = 0;
  for (let i = 0; i < layoutJson.length; i += MAX_STORAGE_SIZE) {
    const part = layoutJson.substring(i, i + MAX_STORAGE_SIZE);
    const key = `${studyId}|${viewName}|${partNumber}`;
    await grok.userSettings.add(STORAGE_NAME, key, part);
    partNumber++;
  }

  if (totalParts > 1) {
    const metadataKey = `${studyId}|${viewName}|metadata`;
    await grok.userSettings.add(STORAGE_NAME, metadataKey, JSON.stringify({parts: totalParts}));
  }
}

export async function loadLayoutFromUserStorage(studyId: string, viewName: string): Promise<string | undefined> {
  try {
    const metadataKey = `${studyId}|${viewName}|metadata`;
    const metadataStr = await grok.userSettings.getValue(STORAGE_NAME, metadataKey);

    if (metadataStr) {
      const metadata = JSON.parse(metadataStr);
      const parts: string[] = [];

      for (let partNumber = 0; partNumber < metadata.parts; partNumber++) {
        const key = `${studyId}|${viewName}|${partNumber}`;
        const part = await grok.userSettings.getValue(STORAGE_NAME, key);
        if (!part)
          return;
        parts.push(part);
      }

      return parts.join('');
    } else {
      const key = `${studyId}|${viewName}|0`;
      return await grok.userSettings.getValue(STORAGE_NAME, key);
    }
  }
  catch (e) {
    return;
  }
}

async function deleteLayoutFromUserStorage(studyId: string, viewName: string): Promise<void> {
  try {
    const metadataKey = `${studyId}|${viewName}|metadata`;
    const metadataStr = await grok.userSettings.getValue(STORAGE_NAME, metadataKey);

    if (metadataStr) {
      const metadata = JSON.parse(metadataStr);
      for (let partNumber = 0; partNumber < metadata.parts; partNumber++) {
        const key = `${studyId}|${viewName}|${partNumber}`;
        await grok.userSettings.delete(STORAGE_NAME, key);
      }
      await grok.userSettings.delete(STORAGE_NAME, metadataKey);
    } else {
      const key = `${studyId}|${viewName}|0`;
      await grok.userSettings.delete(STORAGE_NAME, key);
    }
  }
  catch (e) {
    // Ignore errors when deleting
  }
}

export async function setupTableViewLayout(tableView: DG.TableView, studyId: string, viewName: string): Promise<void> {
  const ribbons = tableView.getRibbonPanels();

  const saveLayoutButton = ui.button('Save layout', async () => {
    try {
      const layout = tableView.saveLayout();
      const layoutJson = layout.toJson();
      await saveLayoutToUserStorage(studyId, viewName, layoutJson);
      grok.shell.info('Layout saved successfully');
    }
    catch (e: any) {
      grok.shell.error(`Failed to save layout: ${e?.message ?? e}`);
    }
  });
  ribbons.push([saveLayoutButton]);
  tableView.setRibbonPanels(ribbons);

  try {
    const savedLayoutJson = await loadLayoutFromUserStorage(studyId, viewName);
    if (savedLayoutJson) {
      const layout = DG.ViewLayout.fromJson(savedLayoutJson);
      tableView.loadLayout(layout);
    }
  }
  catch (e: any) {
    console.log(`No saved layout found or failed to load: ${e?.message ?? e}`);
  }
}
