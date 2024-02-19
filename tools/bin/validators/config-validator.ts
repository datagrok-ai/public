import {ValidationResult} from '../utils/interfaces';


export function validateConf(config: any): ValidationResult {
  const vr: ValidationResult = {
    message: 'Validation Error\n\nRun `grok config` to see the config file or ' +
      '`grok config --reset` to restore the default template.\n\nDetails:\n',
    value: false,
    warnings: [],
  };

  if (typeof config !== 'object' || config === null) {
    vr.message += 'Couldn\'t find server configuration in the file.';
    return vr;
  }

  if (!config.hasOwnProperty('default') || !config.default) {
    vr.message += 'The default server is not specified.';
    return vr;
  }

  if (!config.hasOwnProperty('servers') || !config.servers) {
    vr.message += 'The servers information is missing.';
    return vr;
  }

  const defaultServer = config.default;
  let hasDefault = false;
  if (typeof defaultServer !== 'string') {
    vr.message += 'The "default" field should contain a server alias.';
    return vr;
  }

  for (const server in config.servers) {
    const info = config.servers[server];
    if (typeof info !== 'object' || info === null || !('url' in info) || !('key' in info)) {
      vr.message += `Configuration for server "${server}" has the wrong format. Add properties "url" and "key".`;
      return vr;
    }
    // For cases where the values were skipped during `grok config` setup
    if ('url' in info && !info.url) vr.warnings!.push(`Missing URL for server "${server}".`);
    if ('key' in info && !info.key) vr.warnings!.push(`Missing key for server "${server}".`);
    if (server == defaultServer) hasDefault = true;
  }

  if (hasDefault) {
    vr.value = true;
    vr.message = '';
  } else 
    vr.message += `The default server "${defaultServer}" is not in the server list.`;
  

  return vr;
}
