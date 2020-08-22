#! /usr/bin/env node
const inquirer = require('inquirer');
const fs = require('fs');
const path = require('path');

let dir = path.basename(process.cwd());

const QUESTIONS = [
    {
        name: 'package-name',
        type: 'input',
        message: `Package name:`,
        default: dir,
        validate: function (input) {
            if (/^([A-Za-z\-\_\d])+$/.test(input)) return true;
            else return 'Project name may only include letters, numbers, underscores and hashes.';
        }
    },
    {
        name: 'server',
        type: 'input',
        message: 'Server:',
        default: 'https://dev.datagrok.ai',
        validate: function (input) {
            if (/^https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=.]{1,256}([-a-zA-Z0-9()@:%_\+.~#?&//=]*)/.test(input))
                return true;
            else
                return 'Server must be a valid URL';
        }
    },
    {
        name: 'dev-key',
        type: 'input',
        message: function (answers) {
            return `Developer key (get it from ${answers["server"]}/u):`;
        },
        validate: function (input) {
            if (/^([A-Za-z\d\-])+$/.test(input)) return true;
            else return 'Developer key may only include letters, numbers or hyphens';
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
            contents = contents.replace(/#{PACKAGE_NAME}/g, answers['package-name']);
            contents = contents.replace(/#{PACKAGE_NAME_LOWERCASE}/g, answers['package-name'].toLowerCase());
            contents = contents.replace(/#{REMOTE_URL}/g, answers['server'] + '/api');
            contents = contents.replace(/#{REMOTE_KEY}/g, answers['dev-key']);
            if (file === 'npmignore')
                file = '.npmignore';
            if (file === 'gitignore')
                file = '.gitignore';
            fs.writeFileSync(`${packageDir}/${file}`, contents, 'utf8');
        } else if (stats.isDirectory()) {
            fs.mkdirSync(`${packageDir}/${file}`);

            // recursive call
            createDirectoryContents(answers,`${templateDir}/${file}`, `${packageDir}/${file}`);
        }
    });
}