#!/usr/bin/env node
const { spawn } = require("child_process");

let argv = require('minimist')(process.argv.slice(2));

let mode = argv['_'][1];
let host = argv['_'][0];

if (mode !== 'debug' && mode !=='deploy') {
    console.log('Mode must be either debug or deploy');
    return;
}

let subfolder = `${__dirname}\\win\\`;
if (process.platform === "darwin")
    subfolder = `${__dirname}/macos/`;

const ls = spawn(`${subfolder}dart`, [`${subfolder}upload.dart`, '-p', '.', '-r', host, mode]);

ls.stdout.on("data", data => {
    console.log(`${data}`);
});

ls.stderr.on("data", data => {
    console.log(`stderr: ${data}`);
});

ls.on('error', (error) => {
    console.log(`error: ${error.message}`);
});
