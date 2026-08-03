import fs from 'fs';
import path from 'path';
import {help} from '../utils/ent-helpers';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import {generateDomainClients} from './api';

/** The domain-ui library version a scaffolded domain app depends on. */
const domainUiDependency = '^0.1.0';

export function add(args: { _: string[], domain?: string | boolean }) {
  // `--domain` is the only option any `add` entity takes (`grok add app --domain`).
  const nOptions = Object.keys(args).length - 1 - (args.domain === undefined ? 0 : 1);
  const nArgs = args['_'].length;
  if (nArgs < 2 || nArgs > 5 || nOptions > 0) return false;
  const entity = args['_'][1];
  if (args.domain !== undefined && entity !== 'app')
    return color.error('`--domain` applies to `grok add app` only');

  const curDir = process.cwd();
  const curFolder = path.basename(curDir);
  const srcDir = path.join(curDir, 'src');
  const jsPath = path.join(srcDir, 'package.js');
  const tsPath = path.join(srcDir, 'package.ts');
  const detectorsPath = path.join(curDir, 'detectors.js');
  const webpackConfigPath = path.join(curDir, 'webpack.config.js');
  const scriptsDir = path.join(curDir, 'scripts');
  const queryDir = path.join(curDir, 'queries');
  const queryPath = path.join(queryDir, 'queries.sql');
  const connectDir = path.join(curDir, 'connections');
  const packagePath = path.join(curDir, 'package.json');
  const templateDir = path.join(path.dirname(path.dirname(__dirname)));

  // Package directory check
  if (!fs.existsSync(packagePath)) return color.error('`package.json` not found');
  try {
    JSON.parse(fs.readFileSync(packagePath, {encoding: 'utf-8'}));
  } catch (error) {
    color.error(`Error while reading ${packagePath}:`);
    console.error(error);
  }

  // TypeScript package check
  const ts = fs.existsSync(path.join(curDir, 'tsconfig.json'));
  const packageEntry = ts ? tsPath : jsPath;
  const ext = ts ? '.ts' : '.js';

  function validateName(name: string) {
    if (!/^([A-Za-z])+([A-Za-z\d])*$/.test(name)) 
      return color.error('The name may only include letters and numbers. It cannot start with a digit');
    
    return true;
  }

  function insertName(name: string, data: string) {
    for (const repl of ['NAME', 'NAME_TITLECASE', 'NAME_LOWERCASE']) 
      data = utils.replacers[repl](data, name);
    
    return data;
  }

  function createPackageEntryFile() {
    if (!fs.existsSync(srcDir)) fs.mkdirSync(srcDir);
    if (!fs.existsSync(packageEntry)) {
      const contents = fs.readFileSync(path.join(templateDir,
        'package-template', 'src', 'package.js'), 'utf8');
      fs.writeFileSync(packageEntry, contents, 'utf8');
    }
  }

  /** The manifest at [manifestPath], or null when it is missing or unreadable. */
  function readManifest(manifestPath: string): any {
    if (!fs.existsSync(manifestPath)) return null;
    try {
      const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
      if (typeof manifest?.name !== 'string' || manifest?.tables == null) {
        color.error(`${manifestPath} is not a domain schema manifest (no \`name\` / \`tables\`)`);
        return null;
      }
      return manifest;
    } catch (error) {
      color.error(`Error while reading ${manifestPath}:`);
      console.error(error);
      return null;
    }
  }

  /** Copies a schema manifest into `databases/<schema>/` (where `grok api` and the
   * publisher look for it); returns its schema name and table names. */
  function copyManifest(manifestPath: string): [string, string[]] | null {
    const manifest = readManifest(path.resolve(curDir, manifestPath));
    if (manifest == null) return null;
    const target = path.join(curDir, 'databases', manifest.name, 'schema.json');
    if (path.resolve(curDir, manifestPath) !== target) {
      fs.mkdirSync(path.dirname(target), {recursive: true});
      fs.copyFileSync(path.resolve(curDir, manifestPath), target);
      console.log(`Copied the schema manifest to databases/${manifest.name}/schema.json`);
    }
    return [manifest.name, Object.keys(manifest.tables)];
  }

  /** Adds `import {DomainAppView} from '@datagrok-libraries/domain-ui';` to the package
   * entry file, after its existing imports (which webpack's externals depend on). */
  function addDomainUiImport() {
    const contents = fs.readFileSync(packageEntry, 'utf8');
    if (contents.includes('@datagrok-libraries/domain-ui')) return;
    const eol = contents.includes('\r\n') ? '\r\n' : '\n';
    const line = `import {DomainAppView} from '@datagrok-libraries/domain-ui';`;
    const lines = contents.split(/\r?\n/);
    let last = -1;
    for (let i = 0; i < lines.length; i++)
      if (/^(import .*|export .* from .*);\s*$/.test(lines[i])) last = i;
    lines.splice(last + 1, 0, line);
    fs.writeFileSync(packageEntry, lines.join(eol), 'utf8');
  }

  /** `grok add app [name] --domain <schema>[.<table>] | <schema.json path>` — a working
   * browse/CRUD app per table, from the `domain-ui` defaults alone. */
  function addDomainApp(domainArg: string | boolean, appName: string | null): boolean | void {
    if (typeof domainArg !== 'string' || domainArg === '')
      return color.error('`--domain` needs `<schema>` or `<schema>.<table>`, ' +
        'or a path to a schema.json manifest');

    let schema: string;
    let tables: string[];
    if (/\.json$/i.test(domainArg)) {
      const copied = copyManifest(domainArg);
      if (copied == null) return false;
      [schema, tables] = copied;
    } else {
      if (!/^[A-Za-z_]\w*(\.[A-Za-z_]\w*)?$/.test(domainArg))
        return color.error(`Malformed domain address: ${domainArg}. ` +
          'Use `<schema>`, `<schema>.<table>`, or a path to a schema.json manifest');
      const dot = domainArg.indexOf('.');
      schema = dot === -1 ? domainArg : domainArg.slice(0, dot);
      if (dot !== -1)
        tables = [domainArg.slice(dot + 1)];
      else {
        // Table names of a whole schema are only known offline from its manifest.
        const manifest = readManifest(path.join(curDir, 'databases', schema, 'schema.json'));
        if (manifest == null)
          return color.error(`No \`databases/${schema}/schema.json\` in this package: name a table ` +
            `(\`--domain ${schema}.<table>\`) or pass the path to the schema manifest`);
        tables = Object.keys(manifest.tables);
      }
    }
    if (tables.length === 0)
      return color.error(`Schema \`${schema}\` declares no tables`);
    if (appName != null && tables.length > 1)
      return color.error(`\`--domain ${domainArg}\` covers ${tables.length} tables — ` +
        'name a single table to name its app');

    createPackageEntryFile();
    const template = fs.readFileSync(path.join(templateDir, 'entity-template', 'domain-app' + ext), 'utf8');
    const entry = fs.readFileSync(packageEntry, 'utf8');
    const added: string[] = [];
    const addedTables: string[] = [];
    for (const table of tables) {
      // No 'App' postfix — `grok check` asks apps to drop it (check.ts:350-356).
      const funcName = appName ?? utils.snakeToCamelCase(table);
      if (!validateName(funcName)) return false;
      if (new RegExp(`function\\s+${funcName}\\s*\\(`).test(entry)) {
        color.warn(`${funcName} is already declared in ${path.relative(curDir, packageEntry)} — skipped`);
        continue;
      }
      fs.appendFileSync(packageEntry,
        utils.replacers['DOMAIN_TABLE'](insertName(funcName, template), `${schema}.${table}`));
      added.push(funcName);
      addedTables.push(`${schema}.${table}`);
    }
    if (added.length === 0) return true;

    addDomainUiImport();

    const packageObj = JSON.parse(fs.readFileSync(packagePath, 'utf8'));
    packageObj.dependencies = packageObj.dependencies ?? {};
    if (packageObj.dependencies['@datagrok-libraries/domain-ui'] === undefined) {
      packageObj.dependencies['@datagrok-libraries/domain-ui'] = domainUiDependency;
      fs.writeFileSync(packagePath, JSON.stringify(packageObj, null, 2), 'utf8');
    }

    // Manifests in the package get typed clients AND the typed UI wrappers the app
    // code can switch to; a schema deployed by another package has none, and the
    // reflective app above needs none.
    if (!generateDomainClients(curDir, {ui: true})) return false;

    console.log(help.domainApp(added, addedTables));
    return true;
  }

  let name;
  let tag;
  let contents;

  switch (entity) {
  case 'script':
    if (nArgs < 4 || nArgs > 5) return false;
    let lang = args['_'][2];
    name = args['_'][3];
    if (nArgs === 5) {
      tag = args['_'][2];
      lang = args['_'][3];
      name = args['_'][4];
    }
    if (!Object.keys(utils.scriptLangExtMap).includes(lang)) {
      color.error(`Unsupported language: ${lang}`);
      console.log('You can add a script in one of the following languages:');
      console.log(Object.keys(utils.scriptLangExtMap).join(', '));
      return false;
    }

    // Script name check
    if (!validateName(name)) return false;

    if (tag && tag !== 'panel')
      return color.error('Currently, you can only add the `panel` tag');

    // Create the folder `scripts` if it doesn't exist yet
    if (!fs.existsSync(scriptsDir)) fs.mkdirSync(scriptsDir);

    const scriptPath = path.join(scriptsDir, `${name}.${utils.scriptLangExtMap[lang]}`);
    if (fs.existsSync(scriptPath))
      return color.error(`The file with the script already exists: ${scriptPath}`);

    // Copy the script template
    const templatePath = path.join(templateDir, 'script-template', `${lang}.${utils.scriptLangExtMap[lang]}`);
    contents = fs.readFileSync(templatePath, 'utf8');
    if (tag) {
      const isJs = ['javascript', 'nodejs'].includes(lang);
      const comment = isJs ?
        '//meta.role: panel' :
        '#meta.role: panel';

      const nameLineRe = /^(#|\/\/)\s*name:.*$/m;

      if (nameLineRe.test(contents)) 
        contents = contents.replace(nameLineRe, (match) => `${match}\n${comment}`);
      else 
        contents = `${comment}\n${contents}`;
    }

    fs.writeFileSync(scriptPath, insertName(name, contents), 'utf8');

    // Provide a JS wrapper for the script
    console.log(help.script(name, curFolder));
    break;

  case 'app':
    if (args.domain !== undefined) {
      if (nArgs > 3) return false;
      return addDomainApp(args.domain, nArgs === 3 ? args['_'][2] : null);
    }
    if (nArgs !== 3) return false;

    // App name check
    name = args['_'][2];
    if (!validateName(name)) return false;

    // Create src/package.js if it doesn't exist yet
    createPackageEntryFile();

    // Add an app template to package.js
    const app = fs.readFileSync(path.join(templateDir, 'entity-template', 'app.js'), 'utf8');
    fs.appendFileSync(packageEntry, insertName(name, app));
    console.log(help.app(name));
    break;

  case 'function':
    if (nArgs < 3 || nArgs > 4) return false;

    name = args['_'][2];
    if (nArgs === 4) {
      tag = args['_'][2];
      name = args['_'][3];
    }

    if (!validateName(name)) return false;

    if (tag && tag !== 'panel' && tag !== 'init') 
      return color.error('Currently, you can only add the `panel` or `init` tag');
      

    // Create src/package.js if it doesn't exist yet
    createPackageEntryFile();

    // Add a function to package.js
    const filename = tag === 'panel' ? 'panel' + ext :
      tag === 'init' ? 'init.js' : 'function' + ext;
    const func = fs.readFileSync(path.join(templateDir, 'entity-template', filename), 'utf8');
    fs.appendFileSync(packageEntry, insertName(name, func));

    console.log(help.func(name, tag === 'panel'));
    break;
  case 'connection':
    if (nArgs !== 3) return false;
    name = args['_'][2];

    // Connection name check
    if (!validateName(name)) return false;

    // Create the `connections` folder if it doesn't exist yet
    if (!fs.existsSync(connectDir)) fs.mkdirSync(connectDir);

    const connectPath = path.join(connectDir, `${name}.json`);
    if (fs.existsSync(connectPath)) 
      return color.error(`The connection file already exists: ${connectPath}`);
      

    const connectionTemplate = fs.readFileSync(path.join(templateDir,
      'entity-template', 'connection.json'), 'utf8');
    fs.writeFileSync(connectPath, insertName(name, connectionTemplate), 'utf8');
    console.log(help.connection(name));
    break;

  case 'query':
    if (nArgs !== 3) return false;

    // Query name check
    name = args['_'][2];
    if (!validateName(name)) return false;

    // Create the `queries` folder if it doesn't exist yet
    if (!fs.existsSync(queryDir)) fs.mkdirSync(queryDir);

    // Add a query to queries.sql
    const query = fs.readFileSync(path.join(templateDir, 'entity-template', 'queries.sql'), 'utf8');
    contents = insertName(name, query);
    let connection;
    if (fs.existsSync(connectDir) && fs.readdirSync(connectDir).length !== 0) {
      // Use the name of the first found connection
      connection = fs.readdirSync(connectDir).find((c) => /.+\.json$/.test(c))?.slice(0, -5);
    } else {
      // Create the default connection file
      if (!fs.existsSync(connectDir)) fs.mkdirSync(connectDir);
      connection = fs.readFileSync(path.join(templateDir, 'entity-template', 'connection.json'), 'utf8');
      fs.writeFileSync(path.join(connectDir, 'connection.json'), insertName('connection', connection), 'utf8');
      connection = 'connection';
    }
    contents = contents.replace('#{CONNECTION}', connection!);
    fs.appendFileSync(queryPath, contents);
    console.log(help.query(name));
    break;

  case 'view':
    if (nArgs !== 3) return false;

    // View name check
    name = args['_'][2];
    if (!validateName(name)) return false;
    if (!name.endsWith('View')) 
      color.warn('For consistency reasons, we recommend postfixing classes with \'View\'');
      

    // Create src/package.js if it doesn't exist yet
    createPackageEntryFile();

    // Add a new JS file with a view class
    const viewPath = path.join(srcDir, utils.camelCaseToKebab(name) + ext);
    if (fs.existsSync(viewPath)) 
      return color.error(`The view file already exists: ${viewPath}`);
      
    const viewClass = fs.readFileSync(path.join(templateDir, 'entity-template', 'view-class' + ext), 'utf8');
    fs.writeFileSync(viewPath, insertName(name, viewClass), 'utf8');

    // Add a view function to package.js
    const view = fs.readFileSync(path.join(templateDir, 'entity-template', 'view.js'), 'utf8');
    contents = insertName(name, `import {#{NAME}} from './${utils.camelCaseToKebab(name)}';\n`);
    contents += fs.readFileSync(packageEntry, 'utf8');
    contents += insertName(name, view);
    fs.writeFileSync(packageEntry, contents, 'utf8');
    console.log(help.view(name));
    break;

  case 'viewer':
    if (nArgs !== 3) return false;

    // Viewer name check
    name = args['_'][2];
    if (!validateName(name)) return false;
    if (!name.endsWith('Viewer')) 
      color.warn('For consistency reasons, we recommend postfixing classes with \'Viewer\'');
      

    // Create src/package.js if it doesn't exist yet
    createPackageEntryFile();

    // Add a new JS file with a viewer class
    const viewerPath = path.join(srcDir, utils.camelCaseToKebab(name) + ext);
    if (fs.existsSync(viewerPath)) 
      return color.error(`The viewer file already exists: ${viewerPath}`);
      
    const viewerClass = fs.readFileSync(path.join(templateDir,
      'entity-template', 'viewer-class' + ext), 'utf8');
    fs.writeFileSync(viewerPath, insertName(name, viewerClass), 'utf8');

    console.log(help.viewer(name));
    break;
  case 'detector':
    if (nArgs !== 3) return false;
    name = args['_'][2];

    if (!fs.existsSync(detectorsPath)) {
      let temp = fs.readFileSync(path.join(templateDir, 'package-template', 'detectors.js'), 'utf8');
      // eslint-disable-next-line new-cap
      temp = utils.replacers['PACKAGE_DETECTORS_NAME'](temp, curFolder);
      fs.writeFileSync(detectorsPath, temp, 'utf8');
    }

    const detector = fs.readFileSync(path.join(templateDir, 'entity-template', 'sem-type-detector.js'), 'utf8');
    contents = fs.readFileSync(detectorsPath, 'utf8');
    const idx = contents.search(/(?<=PackageDetectors extends DG.Package\s*{\s*(\r\n|\r|\n)).*/);
    if (idx === -1) return color.error('Detectors class not found'); 
    contents = contents.slice(0, idx) + detector + contents.slice(idx);

    for (const repl of ['NAME', 'NAME_PREFIX', 'PACKAGE_DETECTORS_NAME']) 
      contents = utils.replacers[repl](contents, name);
      

    fs.writeFileSync(detectorsPath, contents, 'utf8');
    console.log(help.detector(name));
    break;
  case 'tests':
    if (!fs.existsSync(tsPath)) {
      color.error('Only TypeScript packages are supported');
      return false;
    }

    if (!fs.existsSync(webpackConfigPath)) {
      color.error('Webpack configuration not found');
      return false;
    }

    const config = fs.readFileSync(webpackConfigPath, 'utf8');
    if (!/(?<=entry:\s*{\s*(\r\n|\r|\n))[^}]*test:/.test(config)) {
      const entryIdx = config.search(/(?<=entry:\s*{\s*(\r\n|\r|\n)).*/);
      if (entryIdx === -1)
        return color.error('Entry point not found during config parsing');
  
      const testEntry = '    test: {filename: \'package-test.js\', library: ' +
          '{type: \'var\', name:`${packageName}_test`}, import: \'./src/package-test.ts\'},';
      fs.writeFileSync(webpackConfigPath, config.slice(0, entryIdx) + testEntry +
          config.slice(entryIdx), 'utf8');
    }

    const packageObj = JSON.parse(fs.readFileSync(packagePath, 'utf8'));
    Object.assign(packageObj.dependencies, {
      '@datagrok-libraries/utils': 'latest',
    });    
    Object.assign(packageObj.scripts, {
      'test': 'grok test',
    });
    fs.writeFileSync(packagePath, JSON.stringify(packageObj, null, 2), 'utf8');

    const packageTestPath = path.join(srcDir, 'package-test.ts');
    if (!fs.existsSync(packageTestPath)) {
      fs.writeFileSync(packageTestPath, fs.readFileSync(path.join(templateDir,
        'package-template', 'src', 'package-test.ts')));
    }

    const testsDir = path.join(srcDir, 'tests');
    if (!fs.existsSync(testsDir))
      fs.mkdirSync(testsDir);

    fs.writeFileSync(path.join(testsDir, 'test-examples.ts'),
      fs.readFileSync(path.join(templateDir, 'entity-template', 'test.ts')));

    fs.writeFileSync(packageTestPath, `import './tests/test-examples';\n` +
        fs.readFileSync(packageTestPath));
    console.log(help.test(testsDir));
    break;
  default:
    return false;
  }
  return true;
}
