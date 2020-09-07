const fs = require('fs')
const inquirer = require('inquirer')
const os = require('os')
const path = require('path')
const yaml = require('js-yaml')

module.exports = {
  config: config
}

const confTemplate = {
  servers: {
    dev: { url: 'https://dev.datagrok.ai', key: '' },
    public: { url: 'https://public.datagrok.ai', key: '' },
    local: { url: 'https://localhost:8082', key: '' }
  },
  default: 'public'
}

const grokDir = path.join(os.homedir(), '.grok')
const confPath = path.join(grokDir, 'config.yaml')

function validateKey (key) {
  if (/^([A-Za-z\d-])+$/.test(key)) {
    return true
  } else {
    return 'Developer key may only include letters, numbers, or hyphens'
  }
}

function validateConf (config) {
  let valid = false
  if (config.hasOwnProperty('servers') && config.servers) {
    valid = true
    for (const server in config.servers) {
      const validServer = config.servers[server].hasOwnProperty('url') && config.servers[server].hasOwnProperty('key')
      valid = valid && validServer
    }
  }
  return valid
}

function generateKeyQ (server) {
  const question = {
    name: server,
    type: 'input',
    message: `Developer key (get it from https://${server}.datagrok.ai/u):`,
    validate: validateKey
  }
  if (server === 'local') {
    question.message = 'Developer key for https://localhost:8082:'
  }
  return question
}

function config (args) {
  const nOptions = Object.keys(args).length - 1
  if (args['_'].length === 1 && (nOptions < 1 || nOptions === 1 && args.reset)) {
    if (!fs.existsSync(grokDir)) {
      fs.mkdirSync(grokDir)
    }
    if (!fs.existsSync(confPath) || args.reset) {
      fs.writeFileSync(confPath, yaml.safeDump(confTemplate))
    }
    const config = yaml.safeLoad(fs.readFileSync(confPath))
    console.log(`Your config file (${confPath}):`)
    console.log(config)
    if (!validateConf(config)) {
      console.log('The file is corrupted. Please run `grok config --reset` to restore the default template')
      return false
    }
    (async () => {
      try {
        const answers = await inquirer.prompt({
          name: 'edit-config',
          type: 'confirm',
          message: 'Do you want to edit it?',
          default: false
        })
        if (answers['edit-config']) {
          for (const server in config.servers) {
            const question = generateKeyQ(server)
            question.default = config['servers'][server]['key']
            const devKey = await inquirer.prompt(question)
            config['servers'][server]['key'] = devKey[server]
          }
          const defaultServer = await inquirer.prompt({
            name: 'default-server',
            type: 'input',
            message: 'Your default server:',
            validate: function (server) {
              if (server in config.servers) {
                return true
              } else {
                return 'Only one of the specified servers may be chosen as default'
              }
            },
            default: config.default
          })
          config.default = defaultServer['default-server']
          fs.writeFileSync(confPath, yaml.safeDump(config))
        }
      } catch (err) {
        console.error('The file is corrupted. Please run `grok config --reset` to restore the default template')
        console.error(err)
        return false
      }
    })()
    return true
  }
}
