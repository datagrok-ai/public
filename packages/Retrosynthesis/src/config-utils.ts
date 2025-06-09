import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BASE_PATH, CONFIGS_PATH} from './const';
import {_package} from './package';
import {updateRetrosynthesisWidget} from './utils';
import {parse} from 'yaml';

export const STORAGE_NAME = 'retrosynthesis';
export const KEY = 'config';
export const DEFAULT_CONFIG_NAME = 'default';
export const TOKEN_PARAM_NAME = 'token';
export const CONFIG_PARAMS = ['expansion', 'filter', 'stock'];
export interface UserRetrosynthesisConfig {
  configName: string,
  expansion: string,
  filter: string,
  stock: string,
}

export interface RetrosynthesisPolicies {
  expansion: string[],
  filter: string[],
  stock: string[],
}

export let currentUserConfig: UserRetrosynthesisConfig =
  {configName: DEFAULT_CONFIG_NAME, expansion: '', filter: '', stock: ''};

export function configIcon(currentMolecule: string, widget: DG.Widget): HTMLElement {
  const settings = ui.icons.settings(async () => {
    const configsInAppData = await getConfigFilesFromAppData();
    const configChoices = [DEFAULT_CONFIG_NAME].concat(configsInAppData);
    await setValidUserConfig();
    const policiesDiv = ui.div();

    const newConfig = Object.assign({}, currentUserConfig);
    const updateAdditionalParams = async (first?: boolean) => {
      ui.empty(policiesDiv);
      const paramsDiv = ui.divV([]);
      if (configFolderInput.value !== DEFAULT_CONFIG_NAME) {
        const configJson = await getCustomConfigJson(configFolderInput.value!);
        if (configJson) {
          //expansion policy
          if (Object.keys(configJson.expansion).length) {
            const expansionChoice = ui.input.choice('expansion', {
              value: first ? currentUserConfig.expansion : '',
              items: [''].concat(Object.keys(configJson.expansion)),
              onValueChanged: () => {
                newConfig.expansion = expansionChoice.value!;
              },
            });
            paramsDiv.append(expansionChoice.root);
          }

          //stock
          if (Object.keys(configJson.stock).length) {
            const stockChoice = ui.input.choice('stock', {
              value: first ? currentUserConfig.stock : '',
              items: [''].concat(Object.keys(configJson.stock)),
              onValueChanged: () => {
                newConfig.stock = stockChoice.value!;
              },
            });
            paramsDiv.append(stockChoice.root);
          }

          //filter policy
          if (Object.keys(configJson.filter).length) {
            const filterChoice = ui.input.choice('filter', {
              value: first ? currentUserConfig.filter : '',
              items: [''].concat(Object.keys(configJson.filter)),
              onValueChanged: () => {
                newConfig.filter = filterChoice.value!;
              },
            });
            paramsDiv.append(filterChoice.root);
          }
        }
      }
      policiesDiv.append(paramsDiv);
    };

    const configFolderInput = ui.input.choice('Config', {
      value: currentUserConfig.configName,
      items: configChoices,
      onValueChanged: async () => {
        currentUserConfig.configName = configFolderInput.value!;
        updateAdditionalParams();
      },
    });

    updateAdditionalParams(true);
    const dlg = ui.dialog('Settings')
      .add(ui.divV([configFolderInput.root, policiesDiv]))
      .onOK(() => {
        if (Object.keys(currentUserConfig).some((key) => (currentUserConfig as any)[key] !== (newConfig as any)[key])) {
          grok.userSettings.add(STORAGE_NAME, KEY, JSON.stringify(newConfig));
          grok.shell.info(`Current config updated`);
          updateRetrosynthesisWidget(currentMolecule, widget);
        }
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

export async function setValidUserConfig(availableConfigs?: string[]): Promise<void> {
  const storedConfig = getStoredUserConfig();
  if (!availableConfigs)
    availableConfigs = await getConfigFilesFromAppData();
  const validConfig = await updateConfigConsideringConfigYml(storedConfig, availableConfigs);
  if (validConfig.error)
    grok.shell.warning(validConfig.error);
  currentUserConfig = validConfig.config;
}

export function getStoredUserConfig(): UserRetrosynthesisConfig {
  const currentConfigStr = grok.userSettings.getValue(STORAGE_NAME, KEY);
  let storedConfig: UserRetrosynthesisConfig;
  if (!currentConfigStr)
    storedConfig = {configName: DEFAULT_CONFIG_NAME, expansion: '', filter: '', stock: ''};
  else {
    try { //for compatibility with previous user storage, which contained only config name
      storedConfig = JSON.parse(currentConfigStr);
    } catch (e) {
      storedConfig = {configName: currentConfigStr, expansion: '', filter: '', stock: ''};
    }
  }
  return storedConfig;
}

export async function updateConfigConsideringConfigYml(storedConfig: UserRetrosynthesisConfig,
  availableConfigs: string[]): Promise<{config: UserRetrosynthesisConfig, error: string}> {
  const updatedConfig: UserRetrosynthesisConfig =
    {configName: DEFAULT_CONFIG_NAME, expansion: '', filter: '', stock: ''};
  if (!availableConfigs.includes(storedConfig.configName)) {
    return {
      config: updatedConfig,
      error: `Config folder ${storedConfig.configName} not found`,
    };
  } else
    updatedConfig.configName = storedConfig.configName;

  try {
    let error = '';
    const configJson = await getCustomConfigJson(storedConfig.configName);
    if (storedConfig.expansion && !configJson?.expansion.includes(storedConfig.expansion))
      error = `Expansion policy ${storedConfig.expansion} not found in config.yml`;
    else
      updatedConfig.expansion = storedConfig.expansion;

    if (storedConfig.stock && !configJson?.stock.includes(storedConfig.stock))
      error = `Stock ${storedConfig.stock} not found in config.yml`;
    else
      updatedConfig.stock = storedConfig.stock;

    if (storedConfig.filter && !configJson?.filter.includes(storedConfig.filter))
      error = `Filter policy ${storedConfig.filter} not found in config.yml`;
    else
      updatedConfig.filter = storedConfig.filter;

    return {
      config: updatedConfig,
      error: error,
    };
  } catch (e: any) {
    return {
      config: updatedConfig,
      error: e?.message ?? e,
    };
  }
}

export async function getCustomConfigJson(configName: string): Promise< RetrosynthesisPolicies| null> {
  let res = null;
  try {
    const fileContent = await grok.dapi.files
      .readAsText(`${BASE_PATH}/${_package.name}/${CONFIGS_PATH}/${configName}/config.yml`);
    res = parse(fileContent);
  } catch (e: any) {
    grok.shell.error(`Error reading config.yml file: ${e?.message ?? e}`);
  }
  return res;
}

