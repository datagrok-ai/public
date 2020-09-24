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

    // Package directory check
    const curDir = process.cwd();
    const curFolder = path.basename(curDir);
    const srcDir = path.join(curDir, 'src');
    const jsPath = path.join(srcDir, 'package.js');
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

    function validateName(name) {
        if (!/^([A-Za-z])+([A-Za-z\d])*$/.test(name)) {
            return console.log('The name may only include letters and numbers. It cannot start with a digit');
        }
        return true;
    }

    function insertName(name, data) {
        data = data.replace(/#{NAME}/g, name)
                   .replace(/#{NAME_TITLECASE}/g,
                    name[0].toUpperCase() + name.slice(1).toLowerCase())
                   .replace(/#{NAME_LOWERCASE}/g, name.toLowerCase())
        return data;
    }

    function createJsFile() {
        if (!fs.existsSync(srcDir)) fs.mkdirSync(srcDir);
        if (!fs.existsSync(jsPath)) {
            var contents = fs.readFileSync(path.join(path.dirname(path.dirname(__dirname)),
                                           'package-template', 'src', 'package.js'), 'utf8');
            fs.writeFileSync(jsPath, contents, 'utf8');
        }
    }

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

            // Script name check
            var name = args['_'][3];
            if (!validateName(name)) return false;

            const scriptsDir = path.join(curDir, 'scripts');
            if (!fs.existsSync(scriptsDir)) fs.mkdirSync(scriptsDir);

            let scriptPath = path.join(scriptsDir, name + '.' + langs[lang]);
            if (fs.existsSync(scriptPath)) {
                return console.log(`The file with the script already exists: ${scriptPath}`);
            }

            // Copy the script template
            let templatePath = path.join(path.dirname(path.dirname(__dirname)), 'script-template')
            templatePath = path.join(templatePath, lang + '.' + langs[lang]);
            var contents = fs.readFileSync(templatePath, 'utf8');
            fs.writeFileSync(scriptPath, insertName(name, contents), 'utf8');

            // Provide a JS wrapper for the script
            console.log(`The script has been created. To call it from a JavaScript file, use:
            
            await grok.functions.call('${curFolder}:${name}', { params })
            
            Read more at https://datagrok.ai/help/compute/scripting
            See examples at https://public.datagrok.ai/scripts,
            https://public.datagrok.ai/js/samples/scripting/scripting`.replace(/ {4}/g, ''));
            break;

        case 'app':
            if (nArgs !== 3) return false;

            // App name check
            var name = args['_'][2];
            if (!validateName(name)) return false;

            // Create src/package.js if it doesn't exist yet
            createJsFile();

            // Add an app template to package.js
            let app = fs.readFileSync(path.join(path.dirname(path.dirname(__dirname)),
                                      'entity-template', 'app.js'), 'utf8');
            fs.appendFileSync(jsPath, insertName(name, app));
            console.log(`The application ${name} has been added successfully`);
            console.log('Read more at https://datagrok.ai/help/develop/develop#applications');
            console.log('See application examples at https://public.datagrok.ai/apps');
            break;

        case 'function':
            if (nArgs !== 3) return false;

            // Function name check
            var name = args['_'][2];
            if (!validateName(name)) return false;

            // Create src/package.js if it doesn't exist yet
            createJsFile();

            // Add a function to package.js
            let func = fs.readFileSync(path.join(path.dirname(path.dirname(__dirname)),
                                      'entity-template', 'function.js'), 'utf8');
            fs.appendFileSync(jsPath, insertName(name, func));
            console.log(`The function ${name} has been added successfully`);
            console.log('Read more at https://datagrok.ai/help/overview/functions/function');
            console.log('See examples at https://public.datagrok.ai/functions');
            break;
        // case 'group':
        //     break;
        // case 'connection':
        //     break;
        // case 'query':
        //     break;
        // case 'view':
        //     break;
        // case 'panel':
        //     break;
        case 'viewer':
            if (nArgs !== 3) return false;

            // Viewer name check
            var name = args['_'][2];
            if (!validateName(name)) return false;

            // Create src/package.js if it doesn't exist yet
            createJsFile();

            // Add a new JS file with a viewer class
            let viewerPath = path.join(srcDir, `${name.toLowerCase()}.js`);
            if (fs.existsSync(viewerPath)) {
                return console.log(`The viewer file already exists: ${viewerPath}`);
            }
            let viewerClass = fs.readFileSync(path.join(path.dirname(path.dirname(__dirname)),
                                              'entity-template', 'viewer-class.js'), 'utf8');
            fs.writeFileSync(viewerPath, insertName(name, viewerClass), 'utf8');


            // Add a viewer function to package.js
            let viewer = fs.readFileSync(path.join(path.dirname(path.dirname(__dirname)),
                                      'entity-template', 'viewer.js'), 'utf8');
            var contents = `import ${name.toLowerCase()}Viewer from './${name.toLowerCase()}.js'\n`
            contents += fs.readFileSync(jsPath, 'utf8');
            contents += insertName(name, viewer);
            fs.writeFileSync(jsPath, contents);
            console.log(`The viewer ${name} has been added successfully`);
            console.log('Read more at https://datagrok.ai/help/develop/js-api#custom-viewers');
            console.log('See examples at https://github.com/datagrok-ai/public/tree/master/packages/Viewers,');
            console.log('https://public.datagrok.ai/js/samples/functions/custom-viewers/viewers');
            break;
        default:
            return false;
    }
    return true;
}
