const fs = require('fs');
const path = require('path');

module.exports = {
  delete: deletePackage
};

const curDir = process.cwd();
const curFolder = path.basename(curDir);

function deletePackage(args) {
  const nOptions = Object.keys(args).length - 1;
  const nArgs = args['_'].length;

  if (nArgs > 2 || nOptions > 0) return false;

  let name = curFolder;
  if (nArgs === 2) name = args['_'][1];

  let packageDir = curDir;
  if (curFolder !== name) packageDir = path.join(packageDir, name);
  if (!fs.existsSync(packageDir)) return console.log(`${packageDir} not found`);

  let packageFile = path.join(packageDir, 'package.json');
  if (!fs.existsSync(packageFile)) return console.log('`package.json` not found');

  try {
    let package = JSON.parse(fs.readFileSync(packageFile));
    if (package.name !== name.toLowerCase()) {
      return console.log('The package name differs from the one in `package.json`');
    }

    fs.rmdirSync(packageDir, {recursive: true});
  } catch (error) {
    console.error(`Error while deleting ${packageDir}:`)
    console.error(error);
  }

  console.log(`${name} is deleted`);
  return true;
}
