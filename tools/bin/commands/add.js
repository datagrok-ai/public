const fs = require('fs');
const path = require('path');

module.exports = {
    add: add
};

function add(args) {
    const nOptions = Object.keys(args).length - 1;
    const nArgs = args['_'].length;
    if (nArgs < 3 || nArgs > 4 || nOptions > 0) return false;
    const entity = args['_'][1];

    switch (entity) {
        case 'script':
            if (nArgs !== 4) return false;
            const lang = args['_'][2];
            const langs = {javascript: 'js', julia: 'jl',
                           node: 'js', octave: 'm', python: 'py', r: 'R'};
            if (!Object.keys(langs).includes(lang)) {
                console.log('Unsupported language');
                console.log('You can add a script in one of the following languages:');
                console.log(Object.keys(langs).join(', '));
                return false;
            }

            // Package directory check
            const curDir = process.cwd();
            const curFolder = path.basename(curDir);
            const packagePath = path.join(curDir, 'package.json');
            if (!fs.existsSync(packagePath)) return console.log('`package.json` not found');
            try {
                const package = JSON.parse(fs.readFileSync(packagePath));
                if (package.fullName !== curFolder) {
                    return console.log('The package name differs from the one in `package.json`');
                }
            } catch (error) {
                console.error(`Error while reading ${packagePath}:`)
                console.error(error);
            }

            // Script name check
            let name = args['_'][3];
            if (!/^([A-Za-z\d])+$/.test(name)) {
                return console.log('A script name may only include letters and numbers');
            }

            const scriptsDir = path.join(curDir, 'scripts');
            if (!fs.existsSync(scriptsDir)) fs.mkdirSync(scriptsDir);

            let scriptPath = path.join(scriptsDir, name + '.' + langs[lang]);
            if (fs.existsSync(scriptPath)) {
                return console.log(`The file with the script already exists: ${scriptPath}`);
            }

            // Copy the script template
            let templatePath = path.join(path.dirname(path.dirname(__dirname)), 'script-template')
            templatePath = path.join(templatePath, lang + '.' + langs[lang]);
            let contents = fs.readFileSync(templatePath, 'utf8');
            contents = contents.replace('Template', name);
            fs.writeFileSync(scriptPath, contents, 'utf8');
            break;
        // case 'app':
        //     break;
        // case 'viewer':
        //     break;
        default:
            return false;
    }
    return true;
}
