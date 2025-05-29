import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BASE_PATH, CONFIGS_PATH} from './const';
import {_package} from './package';

export const STORAGE_NAME = 'retrosynthesis';
export const KEY = 'config';
export const DEFAULT_CONFIG_NAME = 'default';
export const TOKEN_PARAM_NAME = 'token';
let token = '';

export function configIcon(): HTMLElement {
  const settings = ui.icons.settings(async () => {
    const currentConfig = grok.userSettings.getValue(STORAGE_NAME, KEY);
    const configsInAppData = await getConfigFilesFromAppData();
    const configChoices = [DEFAULT_CONFIG_NAME].concat(configsInAppData);
    const configFolderInput = ui.input.choice('Config', {
      value: currentConfig && configChoices.includes(currentConfig) ? currentConfig : DEFAULT_CONFIG_NAME,
      items: configChoices,
    });
    const dlg = ui.dialog('Settings')
      .add(configFolderInput.root)
      .onOK(() => {
        grok.userSettings.add(STORAGE_NAME, KEY,
          configFolderInput!.value! === DEFAULT_CONFIG_NAME ? '' : configFolderInput!.value!);
        grok.shell.info(`Current config saved`);
      });
    dlg.root.classList.add('retrosynthesis-settings-dlg');
    dlg.show();
  });
  return settings;
}

export async function getConfigFilesFromAppData(): Promise<string[]> {
  const targetsFiles: DG.FileInfo[] = await grok.dapi.files
    .list(`${BASE_PATH}/${_package.name}/${CONFIGS_PATH}`, false);
  const configDirs = targetsFiles.filter((file) => file.isDirectory).map((dir) => dir.fileName);
  return configDirs;
}

export async function addUserDefinedConfig(file: DG.FileInfo): Promise<void> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const fileContent = await file.readAsString();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;

  console.log(JSON.stringify({config: fileContent, user_id: userId, config_name: file.name}));
  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/add_user_custom_config', {
      method: 'POST',
      body: JSON.stringify({config: fileContent, user_id: userId, config_name: file.name}),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to add configuration: ${errorText}`);
    }

    grok.shell.info('Configuration added successfully');
  } catch (error) {
    grok.shell.error(`Error adding configuration: ${error}`);
    throw error;
  }
}

export async function getUserConfigsFromDocker(): Promise<string[]> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;

  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, `/get_user_configs`, {
      method: 'POST',
      body: JSON.stringify({user_id: userId}),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to get configuration files: ${errorText}`);
    }

    const configs = await response.json();
    return configs.configs ?? [];
  } catch (error) {
    grok.shell.error(`Error getting configuration files: ${error}`);
    throw error;
  }
}

export async function getConfigsFromDocker(): Promise<string[]> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();

  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, `/get_configs`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to get config folders: ${errorText}`);
    }

    const configs = await response.json();
    return configs.result ?? [];
  } catch (error) {
    grok.shell.error(`Error getting config folders: ${error}`);
    throw error;
  }
}

export async function addConfigToDocker(files: string[], folder: string): Promise<void> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();

  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, `/save_config_files`, {
      method: 'POST',
      body: JSON.stringify({files: files, token: token, folder: folder}),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to load config: ${errorText}`);
    }

    const configs = await response.json();
    grok.shell.info('Configuration added successfully');
  } catch (error) {
    grok.shell.error(`Error loading config folder: ${error}`);
    throw error;
  }
}

async function getToken() {
  if (token === '') {
    const credentials = await _package.getCredentials();
    if (!credentials)
      throw new Error('Token is not set in package credentials');
    if (!credentials.parameters[TOKEN_PARAM_NAME])
      throw new Error('Token is not set in package credentials');
    token = credentials.parameters[TOKEN_PARAM_NAME];
  }
}

function getHostName() {
  const host = window.location.origin;
  if (host.includes('://localhost') || host.includes('://127.0.0.1'))
    return `http://host.docker.internal:8082`;
  return `${host}/api`;
}

export async function syncConfig(dirName: string): Promise<void> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  try {
    await getToken();
    const url = getHostName();
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/sync_dir', {
      method: 'POST',
      body: JSON.stringify({
        from_dir_name: `${_package.name}/${dirName}`,
        to_dir_name: dirName,
        token: token,
        url: url,
      }),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to sync config: ${errorText}`);
    }

    grok.shell.info('Configuration synchronized successfully');
  } catch (error) {
    grok.shell.error(`Failed to sync config: ${error}`);
    throw error;
  }
}
