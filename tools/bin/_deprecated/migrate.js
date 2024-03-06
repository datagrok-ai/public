const fs = require('fs');
const path = require('path');
const os = require('os');
const yaml = require('js-yaml');
const utils = require('../utils/utils.js');

module.exports = {
  migrate: migrate,
};

const curDir = process.cwd();
const keysDir = path.join(curDir, 'upload.keys.json');
const packDir = path.join(curDir, 'package.json');
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');
const confTemplate = yaml.load(fs.readFileSync(confTemplateDir));

const grokMap = {
  'datagrok-upload': 'grok publish',
  'debug': '',
  'deploy': '--release',
  'build': '',
  'rebuild': '--rebuild',
};

const replRegExp = new RegExp(Object.keys(grokMap).join('|'), 'g');

function migrate(args) {
  const nOptions = Object.keys(args).length - 1;
  const nArgs = args['_'].length;

  if (nArgs > 1 || nOptions > 0) return false;

  // Create `config.yaml` if it doesn't exist yet
  if (!fs.existsSync(grokDir)) fs.mkdirSync(grokDir);
  if (!fs.existsSync(confPath)) fs.writeFileSync(confPath, yaml.dump(confTemplate));

  const config = yaml.load(fs.readFileSync(confPath));

  // Copy keys to the `config.yaml` file
  if (fs.existsSync(keysDir)) {
    try {
      const keys = JSON.parse(fs.readFileSync(keysDir));
      const urls = utils.mapURL(config);
      for (const url in keys) {
        try {
          let hostname = (new URL(url)).hostname;
          if (url in urls) hostname = urls[url];
          config['servers'][hostname] = {};
          config['servers'][hostname]['url'] = url;
          config['servers'][hostname]['key'] = keys[url];
        } catch (error) {
          console.log(`Skipping an invalid URL in \`upload.keys.json\`: ${url}`);
        }
      }
      fs.writeFileSync(confPath, yaml.dump(config));
      console.log(`Migrated data from local \`upload.keys.json\` to ${confPath}`);
      fs.unlinkSync(keysDir);
      console.log('Successfully deleted the file');
    } catch (error) {
      console.error(error);
    }
  } else 
    console.log('Unable to locate `upload.keys.json`');
  

  // Rewrite scripts in `package.json`
  if (!fs.existsSync(packDir)) return console.log('`package.json` doesn\'t exist');
  try {
    const _package = JSON.parse(fs.readFileSync(packDir));
    for (const script in _package.scripts) {
      if (!_package['scripts'][script].includes('datagrok-upload')) continue;
      _package['scripts'][script] = _package['scripts'][script].replace(replRegExp, (match) => grokMap[match]);
    }
    fs.writeFileSync(packDir, JSON.stringify(_package, null, '\t'));
    console.log('Converting scripts in `package.json`... Done!');
  } catch (error) {
    console.error(error);
  }

  return true;
}
