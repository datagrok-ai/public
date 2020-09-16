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
const confTemplateDir = path.join(path.dirname(path.dirname(__dirname)), 'config-template.yaml');
const confTemplate = yaml.safeLoad(fs.readFileSync(confTemplateDir));
const curDir = process.cwd();
const packDir = path.join(curDir, 'package.json');

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

function mapURL(conf) {
    let urls = {};
    for (let server in conf.servers) {
        urls[conf['servers'][server]['url']] = server;
    }
    return urls;
}

function publish(args) {
    const nOptions = Object.keys(args).length - 1;
    const nArgs = args['_'].length;

    if (nArgs > 2 || nOptions > 4) return false;
    if (!Object.keys(args).slice(1).every(option =>
        ['build', 'rebuild', 'debug', 'release', 'k', 'key'].includes(option))) return false;
    if ((args.build && args.rebuild) || (args.debug && args.release)) return console.log('You have used incompatible options');
    
    // Create `config.yaml` if it doesn't exist yet
    if (!fs.existsSync(grokDir)) fs.mkdirSync(grokDir);
    if (!fs.existsSync(confPath)) fs.writeFileSync(confPath, yaml.safeDump(confTemplate));

    let config = yaml.safeLoad(fs.readFileSync(confPath));
    let host = config.default;
    let alias = config.default;
    let urls = mapURL(config);
    if (nArgs === 2) host = args['_'][1];

    // The host can be passed either as a URL or an alias
    try {
        host = new URL(host);
        alias = (host in urls) ? urls[host] : host.hostname;
        if (!(alias in config.servers)) {
            config['servers'][alias] = {};
            config['servers'][alias]['url'] = host.href;
            config['servers'][alias]['key'] = args.key || '';
            fs.writeFileSync(confPath, yaml.safeDump(config));
        }
    } catch (error) {
        alias = host;
        if (!(alias in config.servers)) return console.log(`Unknown server alias. Please add it to ${confPath}`);
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
    let url = config['servers'][alias]['url'];

    // Get the package name
    if (!fs.existsSync(packDir)) return console.log('`package.json` doesn\'t exist');
    let package = JSON.parse(fs.readFileSync(packDir));
    let packageName = package.name;

    // Upload the package
    process.env.NODE_TLS_REJECT_UNAUTHORIZED = "0";
    process.on('beforeExit', async () => {
        let code = 0;
        try {
            code = await processPackage(!args.release, Boolean(args.rebuild), url, devKey, packageName)
    
        } catch (error) {
            console.error(error);
            code = 1;
        }
        console.log(`Exiting with code ${code}`);
        process.exit(code);
    });

    return true;
}
