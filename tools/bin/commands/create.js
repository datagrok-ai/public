const fs = require('fs')
const os = require('os')
const path = require('path')
const yaml = require('js-yaml')

module.exports = {
  create: create
}

const curDir = process.cwd()
const curFolder = path.basename(process.cwd())
const grokDir = path.join(os.homedir(), '.grok')
const confPath = path.join(grokDir, 'config.yaml')
const config = yaml.safeLoad(fs.readFileSync(confPath))
const server = config['servers'][config.default]['url']
const key = config['servers'][config.default]['key']

function createDirectoryContents (name, configPath, templateDir, packageDir) {
  const filesToCreate = fs.readdirSync(templateDir)

  filesToCreate.forEach(file => {
    const origFilePath = path.join(templateDir, file)
    const copyFilePath = path.join(packageDir, file)
    const stats = fs.statSync(origFilePath)
    if (stats.isFile()) {
      let contents = fs.readFileSync(origFilePath, 'utf8')
      contents = contents.replace(/#{PACKAGE_NAME}/g, name)
      contents = contents.replace(/#{PACKAGE_NAME_LOWERCASE}/g, name.toLowerCase())
      contents = contents.replace(/#{REMOTE_URL}/g, server + '/api')
      contents = contents.replace(/#{REMOTE_KEY}/g, key)
      if (file === 'npmignore') {
        file = '.npmignore'
      }
      if (file === 'gitignore') {
        file = '.gitignore'
      }
      fs.writeFileSync(copyFilePath, contents, 'utf8')
    } else if (stats.isDirectory()) {
      fs.mkdirSync(copyFilePath)
      // recursive call
      createDirectoryContents(name, configPath, origFilePath, copyFilePath)
    }
  })
}

function create (args) {
  const nOptions = Object.keys(args).length - 1
  if (args['_'].length === 2 && nOptions < 1) {
    const name = args['_'][1]
    const validName = /^([A-Za-z\-_\d])+$/.test(name)
    if (validName) {
      let packageDir = curDir
      if (curFolder !== name && !fs.existsSync(path.join(packageDir, name))) {
        packageDir = path.join(packageDir, name)
        fs.mkdirSync(packageDir)
      }
      const templateDir = path.join(path.dirname(__dirname), '/../', 'package-template')
      createDirectoryContents(name, confPath, templateDir, packageDir)
    } else {
      console.log('Package name may only include letters, numbers, underscores or hyphens')
    }
    return true
  }
}
