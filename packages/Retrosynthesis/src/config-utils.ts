import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const STORAGE_NAME = 'retrosynthesis';
export const KEY = 'config';
export const DEFAULT_CONFIG_NAME = 'default';

export function settingsIcon(configName?: string): HTMLElement {
  let userConfigs: string[] | null = null;
  const settings = ui.icons.settings(async () => {
    const configFileInput = ui.input.file('File', {nullable: true});
    const addConfigIcon = ui.icons.add(async () => {
      if (configFileInput.value) {
        const currentConfig = configChoice.value;
        if (userConfigs!.includes(configFileInput.value.name)) {
          grok.shell.error(`Config ${configFileInput.value.name} already exists`);
          return;
        }
        await addUserDefinedConfig(configFileInput.value);
        userConfigs = userConfigs!.concat([configFileInput.value.name]);
        ui.empty(choicesDiv);
        configChoice = ui.input.choice('Current config', {
          nullable: false,
          value: currentConfig,
          items: [DEFAULT_CONFIG_NAME].concat(userConfigs),
        });
        choicesDiv.append(configChoice.root);
      }
      configFileInput.value = null;
    });
    addConfigIcon.classList.add('retrosynthesis-add-config-icon');
    if (!userConfigs)
      userConfigs = await getUserConfigsFromDocker();
    let configChoice = ui.input.choice('Current config', {
      nullable: false,
      value: configName ?? DEFAULT_CONFIG_NAME,
      items: [DEFAULT_CONFIG_NAME].concat(userConfigs),
    });
    const choicesDiv = ui.div(configChoice.root);
    const settingsDiv = ui.divV([
      choicesDiv,
      ui.divH([configFileInput.root, addConfigIcon]),
    ]);
    const dlg = ui.dialog('Settings')
      .add(settingsDiv)
      .onOK(() => {
        grok.userSettings.add(STORAGE_NAME, KEY,
            configChoice!.value! === DEFAULT_CONFIG_NAME ? '' : configChoice!.value!);
        grok.shell.info(`Current config saved`);
      });
    dlg.root.classList.add('retrosynthesis-settings-dlg');
    dlg.show();
  });
  return settings;
}


export async function addUserDefinedConfig(file: DG.FileInfo): Promise<void> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const fileContent = await file.readAsString();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;

  console.log(JSON.stringify({config: fileContent, user_id: userId, config_name: file.name}));
  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/add_user_config', {
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
