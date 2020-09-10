const fs = require('fs');
const path = require('path');

module.exports = {
    create: create
};

const curDir = process.cwd();
const curFolder = path.basename(curDir);

function createDirectoryContents(name, templateDir, packageDir) {
    const filesToCreate = fs.readdirSync(templateDir);

    filesToCreate.forEach(file => {
        const origFilePath = path.join(templateDir, file);
        const copyFilePath = path.join(packageDir, file);
        const stats = fs.statSync(origFilePath);
        if (stats.isFile()) {
            let contents = fs.readFileSync(origFilePath, 'utf8');
            contents = contents.replace(/#{PACKAGE_NAME}/g, name);
            contents = contents.replace(/#{PACKAGE_NAME_LOWERCASE}/g, name.toLowerCase());
            // In the next version, we do not need the `upload.keys.json` file
            if (file === 'upload.keys.json') {
                return false;
            }
            if (file === 'npmignore') {
                file = '.npmignore';
            }
            if (file === 'gitignore') {
                file = '.gitignore';
            }
            fs.writeFileSync(copyFilePath, contents, 'utf8');
        } else if (stats.isDirectory()) {
            fs.mkdirSync(copyFilePath);
            // recursive call
            createDirectoryContents(name, origFilePath, copyFilePath);
        }
    })
}

function isEmpty(dir) {
    return fs.readdirSync(dir).length === 0;
}

function create(args) {
    const nOptions = Object.keys(args).length - 1;
    const nArgs = args['_'].length;
    if (nArgs <= 2 && nOptions < 1) {
        let name = curFolder;
        if (nArgs === 2) {
            name = args['_'][1];
        }
        const validName = /^([A-Za-z\-_\d])+$/.test(name);
        if (validName) {
            let packageDir = curDir;
            if (curFolder !== name) {
                packageDir = path.join(packageDir, name);
                if (!fs.existsSync(packageDir)) {
                    fs.mkdirSync(packageDir);
                }
            }
            if (!isEmpty(packageDir)) {
                console.log();
                console.log('The package directory should be empty');
                return false;
            }
            const templateDir = path.join(path.dirname(path.dirname(__dirname)), 'package-template');
            createDirectoryContents(name, templateDir, packageDir);
        } else {
            console.log('Package name may only include letters, numbers, underscores, or hyphens');
        }
        return true;
    }
}
