const fs = require('fs');
const path = require('path');
const os = require('os');
const yaml = require('js-yaml');

module.exports = {
    migrate: migrate
};

const curDir = process.cwd();
const keysDir = path.join(curDir, 'upload.keys.json');
const packDir = path.join(curDir, 'package.json');
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');
const confTemplate = yaml.safeLoad(fs.readFileSync(confTemplateDir));

const grokMap = {
    'datagrok-upload': 'grok publish',
    'debug': '',
    'deploy': '--release',
    'build': '',
    'rebuild': '--rebuild'
};

const replRegExp = new RegExp(Object.keys(grokMap).join("|"), "g");

function mapURL(conf) {
    let urls = {};
    for (let server in conf.servers) {
        urls[conf['servers'][server]['url']] = server;
    }
    return urls;
}

function migrate(args) {
    const nOptions = Object.keys(args).length - 1;
    const nArgs = args['_'].length;

    if (nArgs > 1 || nOptions > 0) return false;

    // Create `config.yaml` if it doesn't exist yet
    if (!fs.existsSync(grokDir)) fs.mkdirSync(grokDir);
    if (!fs.existsSync(confPath)) fs.writeFileSync(confPath, yaml.safeDump(confTemplate));

    let config = yaml.safeLoad(fs.readFileSync(confPath));

    // Copy keys to the `config.yaml` file
    if (fs.existsSync(keysDir)) {
        try {
            const keys = JSON.parse(fs.readFileSync(keysDir));
            let urls = mapURL(config);
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
            fs.writeFileSync(confPath, yaml.safeDump(config));
            console.log('Migrated data from local `upload.keys.json`');
            fs.unlinkSync(keysDir);
            console.log('Successfully deleted the file');
        } catch (error) {
            console.error(error);
        }
    } else {
        console.log('Unable to locate `upload.keys.json`');
    }

    // Rewrite scripts in `package.json`
    if (!fs.existsSync(packDir)) return console.log('`package.json` doesn\'t exist');
    try {
        let package = JSON.parse(fs.readFileSync(packDir));
        for (let script in package.scripts) {
            if (!package['scripts'][script].includes('datagrok-upload')) continue;
            package['scripts'][script] = package['scripts'][script].replace(replRegExp, (match) => grokMap[match]);
        }
        fs.writeFileSync(packDir, JSON.stringify(package, null, '\t'));
        console.log('Converting scripts in `package.json`... Done!');
    } catch (error) {
        console.error(error);
    }

    return true;
}
