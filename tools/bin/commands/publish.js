const archiver = require('archiver-promise');
const fs = require('fs');
const fetch = require('node-fetch');
const os = require('os');
const path = require('path');
const walk = require('ignore-walk');
const yaml = require('js-yaml');

module.exports = {
    publish: publish
};

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

const curDir = process.cwd();
const keysDir = path.join(curDir, 'upload.keys.json');
const packDir = path.join(curDir, 'package.json');

const grokMap = {
    'datagrok-upload': 'grok publish',
    'debug': '--debug',
    'deploy': '--release',
    'build': '--build',
    'rebuild': '--rebuild'
};

const replRegExp = new RegExp(Object.keys(grokMap).join("|"), "g");

async function processPackage(debug, rebuild, host, devKey, packageName) {
    // Get the server timestamps
    let timestamps = {};
    if (debug) {
        try {
            timestamps = await (await fetch(`${host}/packages/dev/${devKey}/${packageName}/timestamps`)).json();
            if (timestamps['#type'] === 'ApiError') {
                console.log(timestamps.message);
                return 1;
            }
        } catch (error) {
            console.error(error);
            return 1;
        }
    }

    let zip = archiver('zip', {store: false});
    
    // Gather the files
    let localTimestamps = {};
    let files = await walk({
        path: '.',
        ignoreFiles: ['.npmignore', '.gitignore'],
        includeEmpty: true,
        follow: true
    });

    if (!rebuild) {
        const distFiles = await walk({
            path: './dist',
            ignoreFiles: [],
            includeEmpty: true,
            follow: true
        });
        distFiles.forEach((df) => {
            files.push(`dist/${df}`);
        })
    }

    files.forEach((file) => {
        let fullPath = file;
        let relativePath = path.relative(curDir, fullPath);
        let canonicalRelativePath = relativePath.replace(/\\/g, '/');
        if (canonicalRelativePath.includes('/.'))
            return;
        if (canonicalRelativePath.startsWith('.'))
            return;
        if (relativePath.startsWith('node_modules'))
            return;
        if (relativePath.startsWith('dist') && rebuild)
            return;
        if (relativePath.startsWith('upload.keys.json'))
            return;
        if (relativePath === 'zip')
            return;
        let t = fs.statSync(fullPath).mtime.toUTCString();
        localTimestamps[canonicalRelativePath] = t;
        if (debug && timestamps[canonicalRelativePath] === t) {
            console.log(`Skipping ${canonicalRelativePath}`);
            return;
        }
        zip.append(fs.createReadStream(fullPath), {name: relativePath});
        console.log(`Adding ${relativePath}...`);
    });
    zip.append(JSON.stringify(localTimestamps), {name: 'timestamps.json'});
    
    // Upload
    let uploadPromise = new Promise((resolve, reject) => {
            fetch(`${host}/packages/dev/${devKey}/${packageName}?debug=${debug.toString()}&rebuild=${rebuild.toString()}`, {
                method: 'POST',
                body: zip
            }).then(body => body.json()).then(j => resolve(j)).catch(err => {
                reject(err);
            });
        }
    )
    await zip.finalize();

    try {
        let log = await uploadPromise;

        fs.unlinkSync('zip');
        if (log['#type'] === 'ApiError') {
            console.log(log['message']);
            console.log(log['innerMessage']);
            return 1;
        } else {
            console.log(log);
        }
    } catch (error) {
        console.error(error);
        return 1;
    }
    return 0;
}

function publish(args) {
    const nOptions = Object.keys(args).length - 1;
    const nArgs = args['_'].length;

    if (nArgs > 2 || nOptions > 4) return false;
    if (!Object.keys(args).slice(1).every(option =>
        ['build', 'rebuild', 'debug', 'release', 'key', 'migrate'].includes(option))) return false;
    if ((args.build && args.rebuild) || (args.debug && args.release)) return console.log('You have used incompatible options');
    
    if (!fs.existsSync(confPath)) return console.log('There is no `config.yaml` file. Please run `grok config` to create one');
    let config = yaml.safeLoad(fs.readFileSync(confPath));
    
    if (args.migrate) {
        // Copy keys to the `config.yaml` file
        if (fs.existsSync(keysDir)) {
            try {
                const keys = JSON.parse(fs.readFileSync(keysDir));
                for (const url in keys) {
                    try {
                        const hostname = (new URL(url)).hostname;
                        config['servers'][hostname] = {};
                        config['servers'][hostname]['url'] = url.endsWith('/api') ? url.slice(0, -4) : url;
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
        } catch (error) {
            console.error(error);
        }
    }

    config = yaml.safeLoad(fs.readFileSync(confPath));
    let host = config.default;
    if (nArgs === 2) host = args['_'][1];

    // The host can be passed either as a URL or an alias
    try {
        host = new URL(host);
        let alias = host.hostname;
        if (!alias in config.servers) {
            config['servers'][alias] = {};
            config['servers'][alias]['url'] = host.endsWith('/api') ? host.slice(0, -4) : host;
            config['servers'][alias]['key'] = args.key || '';
            fs.writeFileSync(confPath, yaml.safeDump(config));
        }
    } catch (error) {
        alias = host;
        if (!alias in config.servers) return console.log('Unknown server alias. Please add it by running `grok config`');
    }

    // Update the developer key
    if (args.key) {
        config = yaml.safeLoad(fs.readFileSync(confPath));
        config['servers'][alias]['key'] = args.key;
        fs.writeFileSync(confPath, yaml.safeDump(config));
    }

    // Get the non-empty developer key
    config = yaml.safeLoad(fs.readFileSync(confPath));
    let devKey = config['servers'][alias]['key'];
    if (devKey === '') return console.log('Please provide the key with `--key` option or add it by running `grok config`');

    // Get the URL
    let url = new URL(config['servers'][alias]['url']);
    if (url.hostname !== 'localhost' && url.pathname !== '/api') url.href += '/api';

    // Set the modes
    let debug = true;
    if (args.release) debug = false;
    let rebuild = false;
    if (args.rebuild) rebuild = true;

    // Get the package name
    if (!fs.existsSync(packDir)) return console.log('`package.json` doesn\'t exist');
    let package = JSON.parse(fs.readFileSync(packDir));
    let packageName = package.name;

    // Upload the package
    process.env.NODE_TLS_REJECT_UNAUTHORIZED = "0";
    process.on('beforeExit', async () => {
        let code = 0;
        try {
            code = await processPackage(debug, rebuild, url, devKey, packageName)
    
        } catch (error) {
            console.error(error);
            code = 1;
        }
        console.log(`Exiting with code ${code}`);
        process.exit(code);
    });

    return true;
}
