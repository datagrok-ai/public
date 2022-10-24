import fs from 'fs';
import inquirer from 'inquirer';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import { validateConf } from '../validators/config-validator';
import { Config, Indexable } from '../utils/utils';
import * as color from '../utils/color-utils';


const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');
const confTemplate = yaml.load(fs.readFileSync(confTemplateDir, { encoding: 'utf-8' }));

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

function validateKey(key: string) {
  if (!key || /^([A-Za-z\d-])+$/.test(key)) {
    return true;
  } else {
    return 'Developer key may only include letters, numbers, or hyphens';
  }
}

function generateKeyQ(server: string, url: string): Indexable {
  const origin = (new URL(url)).origin;
  const question = {
    name: server,
    type: 'input',
    message: `Developer key (get it from ${origin}/u or press ENTER to skip):`,
    validate: validateKey
  };
  if (server.startsWith('local')) {
    question.message = `Developer key for ${origin} (press ENTER to skip):`;
  }
  return question;
}

async function addNewServer(config: Config) {
  while (true) {
    const addServer = (await inquirer.prompt({
      name: 'add-server',
      type: 'confirm',
      message: 'Do you want to add a server?',
      default: false,
    }))['add-server'];
    if (addServer) {
      const name = (await inquirer.prompt({
        name: 'server-name',
        type: 'input',
        message: 'Enter a name:',
      }))['server-name'];

      const url = (await inquirer.prompt({
        name: 'server-url',
        type: 'input',
        message: 'Enter a URL:',
      }))['server-url'];

      const key = (await inquirer.prompt(generateKeyQ(name, url)) as Indexable)[name];
      config.servers[name] = { url, key };
    } else {
      break;
    }
  }
}

export function config(args: ConfigArgs) {
  const nOptions = Object.keys(args).length - 1;
  const interactiveMode = args['_'].length === 1 && (nOptions < 1 || nOptions === 1 && args.reset);
  const hasAddServerCommand = args['_'].length === 2 && args['_'][1] === 'add' &&
    args.server && args.key && args.k && args.alias && (nOptions === 4 || nOptions === 5 && args.default);
  if (!interactiveMode && !hasAddServerCommand) return false;

  if (!fs.existsSync(grokDir)) {
    fs.mkdirSync(grokDir);
  }
  if (!fs.existsSync(confPath) || args.reset) {
    fs.writeFileSync(confPath, yaml.dump(confTemplate));
  }
  const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as Config;

  if (hasAddServerCommand) {
    try {
      new URL(args.server!);
    } catch (error) {
      color.error('URL parsing error. Please, provide a valid server URL.');
      return false;
    }
    config.servers[args.alias!] = { url: args.server!, key: args.key! };
    color.success(`Successfully added the server to ${confPath}.`);
    console.log(`Use this command to deploy packages: grok publish ${args.alias!}`);
    if (args.default) {
      config.default = args.alias!;
    }
    fs.writeFileSync(confPath, yaml.dump(config));
    return true;
  }

  console.log(`Your config file (${confPath}):`);
  console.log(config);
  const valRes = validateConf(config);
  if (!config || !valRes.value) {
    color.error(valRes.message);
    return false;
  }
  if (valRes.warnings!.length) color.warn(valRes.warnings!.join('\n'));
  (async () => {
    try {
      const answers = await inquirer.prompt({
        name: 'edit-config',
        type: 'confirm',
        message: 'Do you want to edit it?',
        default: false
      });
      if (answers['edit-config']) {
        for (const server in config.servers) {
          const url = config['servers'][server]['url'];
          const question = generateKeyQ(server, url);
          question.default = config['servers'][server]['key'];
          const devKey: Indexable = await inquirer.prompt(question);
          config['servers'][server]['key'] = devKey[server];
        }
        await addNewServer(config);
        const defaultServer = await inquirer.prompt({
          name: 'default-server',
          type: 'input',
          message: 'Your default server:',
          validate: function (server) {
            if (server in config.servers) {
              return true;
            } else {
              return 'Only one of the specified servers may be chosen as default';
            }
          },
          default: config.default
        });
        config.default = defaultServer['default-server'];
        fs.writeFileSync(confPath, yaml.dump(config));
      }
    } catch (err) {
      color.error('The file is corrupted. Please run `grok config --reset` to restore the default template');
      console.error(err);
      return false;
    }
  })();
  fs.writeFileSync(confPath, yaml.dump(config));
  return true;
}

interface ConfigArgs {
  _: string[],
  alias?: string,
  default?: boolean,
  reset?: boolean,
  server?: string,
  key?: string,
  k?: string,
}
