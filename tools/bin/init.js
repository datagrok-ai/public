#! /usr/bin/env node
const inquirer = require('inquirer');
const fs = require('fs');


const QUESTIONS = [
    {
        name: 'package-name',
        type: 'input',
        message: 'Package name:',
        validate: function (input) {
            if (/^([A-Za-z\-\_\d])+$/.test(input)) return true;
            else return 'Project name may only include letters, numbers, underscores and hashes.';
        }
    },
    {
        name: 'remote-server',
        type: 'input',
        message: 'Datagrok server endpoint for debugging (should end with "/api"):',
        validate: function (input) {
            if (/^https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=.]{1,256}([-a-zA-Z0-9()@:%_\+.~#?&//=]*)/.test(input)) return true;
            else return 'Server must be a valid URL';
        }
    },
    {
        name: 'dev-key',
        type: 'input',
        message: 'Datagrok server developer key:',
        validate: function (input) {
            if (/^([A-Za-z\d\-])+$/.test(input)) return true;
            else return 'Datagrok server developer key may only include letters, numbers or hyphens';
        }
    },

];

const CURR_DIR = process.cwd();

inquirer.prompt(QUESTIONS)
    .then(answers => {
        console.log();
        createDirectoryContents(answers, `${__dirname}/../template`, CURR_DIR);
    });

function createDirectoryContents (answers, templateDir, packageDir) {
    const filesToCreate = fs.readdirSync(templateDir);

    filesToCreate.forEach(file => {
        const origFilePath = `${templateDir}/${file}`;

        // get stats about the current file
        const stats = fs.statSync(origFilePath);

        if (stats.isFile()) {
            let contents = fs.readFileSync(origFilePath, 'utf8');
            contents = contents.replace(/#{PACKAGE_NAME}/g, answers['package-name'])
            contents = contents.replace(/#{PACKAGE_NAME_LOWERCASE}/g, answers['package-name'].toLowerCase())
            contents = contents.replace(/#{REMOTE_URL}/g, answers['remote-server'])
            contents = contents.replace(/#{REMOTE_KEY}/g, answers['dev-key'])
            fs.writeFileSync(`${packageDir}/${file}`, contents, 'utf8');
        } else if (stats.isDirectory()) {
            fs.mkdirSync(`${packageDir}/${file}`);

            // recursive call
            createDirectoryContents(answers,`${templateDir}/${file}`, `${packageDir}/${file}`);
        }
    });
}