import { exec } from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import puppeteer from 'puppeteer';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import * as testUtils from '../utils/test-utils';


export function test(args: TestArgs): boolean {
  const options = Object.keys(args).slice(1);
  const commandOptions = ['host', 'csv', 'gui', 'skip-build', 'skip-publish'];
  const nArgs = args['_'].length;
  const curDir = process.cwd();
  const grokDir = path.join(os.homedir(), '.grok');
  const confPath = path.join(grokDir, 'config.yaml');

  if (nArgs > 1 || options.length > commandOptions.length || (options.length > 0 &&
    !options.every(op => commandOptions.includes(op))))
    return false;

  if (!utils.isPackageDir(curDir)) {
    color.error('File `package.json` not found. Run the command from the package directory');
    return false;
  }

  if (!fs.existsSync(confPath)) {
    color.error(`File \`${confPath}\` not found. Run \`grok config\` to set up the config file`);
    return false;
  }

  const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;

  if (args.host) {
    if (args.host in config.servers) {
      process.env.HOST = args.host;
      console.log('Environment variable `HOST` is set to', args.host);
    } else {
      color.error(`Unknown server alias. Please add it to ${confPath}`);
      return false;
    }
  } else if (config.default) {
    process.env.HOST = config.default;
    console.log('Environment variable `HOST` is set to', config.default);
  }

  const packageData = JSON.parse(fs.readFileSync(path.join(curDir, 'package.json'), { encoding: 'utf-8' }));
  if (packageData.name) {
    process.env.TARGET_PACKAGE = utils.kebabToCamelCase(utils.removeScope(packageData.name));
    console.log('Environment variable `TARGET_PACKAGE` is set to', process.env.TARGET_PACKAGE);
  } else {
    color.error('Invalid package name. Set the `name` field in `package.json`');
    return false;
  }

  if (args['skip-build']) {
    if (args['skip-publish'])
      test();
    else
      publish(test);
  } else {
    build(args['skip-publish'] ? test : () => publish(test));
  }

  function build(callback: Function): void {
    exec('npm run build', (err, stdout, stderr) => {
      color.info(`Building package...`);
      if (err) {
        console.log(stdout);
        throw err;
      } else {
        console.log(stdout);
        color.warn(stderr);
      }
      callback();
    });
  }

  function publish(callback: Function): void {
    exec(`grok publish ${process.platform === 'win32' ? '%HOST%' : '${HOST}'}`, (err, stdout, stderr) => {
      color.info(`Publishing package "${process.env.TARGET_PACKAGE}" to ${process.env.HOST}...`);
      if (err)
        throw err;
      else {
        console.log(stdout);
        color.warn(stderr);
      }
      callback();
    });
  }

  function test(): void {
    color.info('Starting tests...');
    const P_START_TIMEOUT: number = 3600000;
    let browser: puppeteer.Browser;
    let page: puppeteer.Page;
    type resultObject = { failReport: string, skipReport: string, passReport: string, failed: boolean, csv?: string };

    function init(timeout: number): Promise<string> {
      const params = Object.assign({}, testUtils.defaultLaunchParameters);
      if (args.gui)
        params['headless'] = false;
      return testUtils.runWithTimeout(timeout, async () => {
        let out = await testUtils.getBrowserPage(puppeteer, params);
        browser = out.browser;
        page = out.page;
        return 'Initialization completed.';
      })
    }
    
    function runTest(timeout: number): Promise<resultObject> {
      return testUtils.runWithTimeout(timeout, async () => {
        const targetPackage: string = process.env.TARGET_PACKAGE ?? '#{PACKAGE_NAMESPACE}';
        console.log(`Testing ${targetPackage} package...\n`);
    
        let r: resultObject = await page.evaluate((targetPackage): Promise<resultObject> => {
          return new Promise<resultObject>((resolve, reject) => {
            (<any>window).grok.functions.eval(targetPackage + ':test()').then((df: any) => {
              let failed = false;
              let skipReport = '';
              let passReport = '';
              let failReport = '';
      
              if (df == null) {
                failed = true;
                failReport = 'Fail reason: No package tests found';
                resolve({failReport, skipReport, passReport, failed});
              }
      
              const cStatus = df.columns.byName('success');
              const cSkipped = df.columns.byName('skipped');
              const cMessage = df.columns.byName('result');
              const cCat = df.columns.byName('category');
              const cName = df.columns.byName('name');
              const cTime = df.columns.byName('ms');
              for (let i = 0; i < df.rowCount; i++) {
                if (cStatus.get(i)) {
                  if (cSkipped.get(i)) {
                    skipReport += `Test result : Skipped : ${cTime.get(i)} : ${targetPackage}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
                  } else {
                    passReport += `Test result : Success : ${cTime.get(i)} : ${targetPackage}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
                  }
                } else {
                  failed = true;
                  failReport += `Test result : Failed : ${cTime.get(i)} : ${targetPackage}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
                }
              }
              const csv = df.toCsv();
              resolve({failReport, skipReport, passReport, failed, csv});
            }).catch((e: any) => reject(e));
          });
        }, targetPackage);
    
        return r;
      });
    }

    (async () => {
      const initMessage = await init(P_START_TIMEOUT);
      console.log(initMessage);
    
      const r = await runTest(7200000);

      if (r.csv && args.csv) {
        fs.writeFileSync(path.join(curDir, 'test-report.csv'), r.csv, 'utf8');
        color.info('Saved `test-report.csv`\n');
      }

      if (r.passReport)
        console.log(r.passReport);

      if (r.skipReport)
        console.log(r.skipReport);

      if (r.failed) {
        console.log(r.failReport);
        color.fail('Tests failed.');
        testUtils.exitWithCode(1);
      } else {
        color.success('Tests passed.');
        testUtils.exitWithCode(0);
      }
    
      //@ts-ignore
      if (browser != null)
        await browser.close();
    })();
  }

  return true;
}

interface TestArgs {
  _: string[],
  host?: string,
  csv?: boolean,
  gui?: boolean,
  'skip-build'?: boolean,
  'skip-publish'?: boolean,
}
