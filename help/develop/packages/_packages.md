<!-- TITLE: Packages -->
<!-- ORDER: 1 -->

# Packages

A Datagrok package is a JavaScript or TypeScript application that introduces a new functionality to Datagrok or extends
the existing one. You can extend Datagrok with _packages_ that you build using the [JavaScript API].

For example, you package might be a [custom viewer built with D3.js], [script for Cheminformatics], and similar. Check
out the [public repository] for more examples of Datagrok packages.

## Table of contents

* [Introduction](#package-elements)
* [Creating a package](#creating-a-package)
* [Package structure](#package-structure)
* [Structuring your package](#structuring-your-package)

## Package elements

A package typically contains the following elements:

* [Functions], [viewers], [widgets], and [applications]
* [Scripts] that can be written in R, Python, Octave, Grok, Julia, JavaScript, Node.js, or Java
* [Database queries] and [connections]

## Creating a package

To develop a package, install the [LTS version of Node.js]. We recommend that you install Node.js using [NVM] to avoid
permission issues when running npm packages globally. Check out the [NPM instructions] for more details.

Now, let's create a package:

1. Install `datagrok-tools` globally to manage your Datagrok packages:

  ```shell
  npm i datagrok-tools -g
  ```

> **Note**: On macOS and Linux, you may also need to install `datagrok-tools` globally with `sudo`: `sudo npm i datagrok-tools -g` .
> When prompted, enter your password.

2. Go to your [user profile on Datagrok] and obtain the developer key.

> **Note**: Datagrok has the [public] and [development servers]. Use the key from the development server. If you're
> using Docker to run a package, refer to the [Deployment with Docker Compose section].

3. Configure your environment:

  ```shell
  grok config
  ```

You will be prompted to enter the developer keys and set the default server. Skip the public key and Docker (local) key.
Datagrok Tools also outputs the `config.yaml` file and the directory it's stored in.

> Your developer keys will be stored locally in `config.yaml` under the home directory &mdash;
> `C:\Users\%Username%\.grok\config.yaml` for Windows and `~/.grok/config.yaml` for Linux. For more details about the
> path to `config.yaml`, refer to the [Node.js > OS Homedir section].

4. Run:

  ```shell
  grok create <package-name>
  ```

> **Note**: If you want to create a TypeScript package, pass the `--ts` option to the command.

5. Cd into the package directory and install other dependencies:

  ```shell
  npm install
  ```

6. After you complete the work on your package, upload it to Datagrok by running the following command:

  ```shell
  grok publish
  ```

Note that by running this command, only you will see the package in Datagrok. For more information about the available
commands, refer to the [Datagrok Tools section].

## Package structure

| File                                    | Description                                               |
|-----------------------------------------|-----------------------------------------------------------|
| node_modules                            | Node.js modules                                           |
| [src/package.js](#package.js)           | The Datagrok package's entry point                        |
| .gitignore                              | Files and folders ignored by Git                          |
| [detectors.js](#detectors)              | File with detectors                                       |
| [package.json](#package.json)           | Package metadata with predefined dependencies and scripts |
| package.PNG                             | Package icon                                              |
| README.md                               | Summary of the package                                    |
| [webpack.config.js](#webpack.config.js) | webpack configuration for development                     |

We discuss the key files in detail in their own sections that follow.

### <a href="#" id="package.json"></a>package.json

The `package.json` file contains the project metadata, including the scripts to build or publish the package,
dependencies, and other data.

```json
{
  "name": "my-grok-app",
  "beta": true,
  "friendlyName": "my-grok-app",
  "version": "0.0.1",
  "description": "",
  "dependencies": {
    "datagrok-api": "latest",
    "cash-dom": "latest",
    "dayjs": "latest"
  },
  "devDependencies": {
    "webpack": "latest",
    "webpack-cli": "latest"
  },
  "scripts": {
    "install-dependencies": "npm install",
    "debug-my-grok-app": "grok publish --rebuild",
    "release-my-grok-app": "grok publish --rebuild --release",
    "build-my-grok-app": "webpack",
    "build": "webpack",
    "debug-my-grok-app-dev": "grok publish dev --rebuild",
    "release-my-grok-app-dev": "grok publish dev --rebuild --release",
    "debug-my-grok-app-local": "grok publish local --rebuild",
    "release-my-grok-app-local": "grok publish local --rebuild --release"
  }
}
```

The package includes the `datagrok-api` dependency &mdash; the JavaScript API you will use to develop the package.
`dayjs` is an [alternative to Moment.js], and `cash-dom` is a [competitor for jQuery]. Additionally, the CLI includes
`webpack` and `webpack-cli` as development dependencies to build the package.

`package.json` also contains a few valuable scripts to let you [debug] or [publish the package]. To run any script in
the terminal, run `npm run <script-name>`.

To install the dependencies, run `npm install` from the terminal or run the provided script `install-dependencies` from
your IDE. You can install other npm packages (such as React or Babel) using `npm install <npm package>`.

### <a href="#" id="package.js"></a>package.js

The `src/package.js` file (or `src/package.ts` if you created a TypeScript package) imports three key Datagrok modules
&mdash; `datagrok-api/grok`, `datagrok-api/ui`, and `datagrok-api/dg`.

```js
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: test
export function test() {
  grok.shell.info(_package.webRoot);
}
```

`package.js` typically contains the [annotated functions] that implement a functionality. Each function needs to be
annotated with respective [header parameters] that define the function metadata. The parameters help Datagrok to
register viewers, widgets, renderers, and other [function types].

If a function has more than one output, it must return a JavaScript object `{prop1: value1, prop2: value2}`, for
example:

```js
//name: test1
//output: string s1
export function test1() {
    // Your logic goes here...
    return 'a';
}
//name: test2
//output: string s1
//output: string s2
export function test2() {
    // Your logic goes here...
    return {s1: 'a', s2: 'b'};
}
```

### <a href="#" id="detectors.js"></a>detectors.js

The `detectors.js` file must define a `<packageName>PackageDetectors` class that extends the base `DG.Package` class.
(For more information about this class, refer to the [Package section] in JavaScript API documentation.)

`detectors.js` is similar to `package.js` but is intended for smaller functions called _[semantic type detectors]_,
which are loaded separately from the rest of the package. Semantic type detectors are called each time a user opens a
table in Datagrok and are used to quickly inspect the data and determine the semantic type of the columns.

Tagging the semantic types allows Datagrok to offer specific functions for data of a particular type.

Here's an example from the [Sequence package] that contains the `detectNucleotides` detector:

```js
class SequencePackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.name.startsWith('nuc')) {
      col.semType = 'nucleotides';
      return 'nucleotides';
    }
    return null;
  }
}
```

Once registered by Datagrok, this function becomes available across the whole platform and can be used for semantic type
detection.

### <a href="#" id="webpack.config.js"></a>webpack.config.js

The Grok CLI generates a typical webpack configuration for the Datagrok package:

```javascript
const path = require('path');

module.exports = {
  // Set 'mode' to 'production' to minify the output and optimize the production build
  mode: 'development',
  entry: {
    // The package is limited to exactly one entry point
    package: './src/package.js'
  },
  // The 'devtool' option enhances package debugging in browser devtools
  devtool: 'inline-source-map',
  // The external modules won't be loaded to the output, but will be taken from the environment
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'dayjs': 'dayjs',
  },
  output: {
    filename: '[name].js',
    // Name of the package in lower case letters
    library: 'sequence',
    // Results will be assigned to a variable `sequence`
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
```

Refer to the [webpack documentation] to modify or extend the provided configuration.

> **Note**:

## Naming conventions

Follow these naming conventions when working on your packages:

* Use upper camel case for package names: `ApiSamples` and `OctaveScripts`
  > **Note**: Package names that comply with the [NPM package naming rules], e.g., `api-samples` and `octave-scripts`,
  > are accepted as well. If you write the desired name in the `friendlyName` field of `package.json`, it will be shown
  > in the UI as the package name.
* Postfix the [views] and [viewers] classes with `~View` and `~Viewer` respectively: `NotebookView`
  , `DumbbellViewer`
* Prefix the names of semantic type detectors with `detect~`: `detectNucleotides`, `detectRDSmiles`
* Use lowercase letters for files, separate words by dashes: `tika-extractor.py`, `chord-viewer.js`

## Structuring your package

Your package might contain the following additional folders, depending on your needs:

| Folder       | Description                            | Documentation               | Package Examples                        |
|--------------|----------------------------------------|-----------------------------|-----------------------------------------|
| environments | Environment configurations             | [Scripting > Environments]  | [DemoScripts]                           |
| scripts      | Collection of scripts for computations | [Scripting]                 | [Chem/scripts], [DemoScripts], [Impute] |
| swaggers     | OpenAPI specifications                 | [Access > OpenAPI]          | [EnamineStore], [Swaggers]              |
| connections  | Database connections                   | [Access > Data Connections] | [Chembl], [UsageAnalysis]               |
| queries      | Database queries                       | [Access > Data Queries]     | [Chembl], [UsageAnalysis]               |
| CSS          | CSS files with custom styles           | -                           | [Notebooks], [Discovery]                |
| data-samples | Data for demo and testing              | -                           | [Chem], [Sunburst]                      |

## What's next?

* [Publishing a package]

[a few valuable scripts]: https://github.com/datagrok-ai/public/tree/master/tools#commands

[alternative to Moment.js]: https://www.npmjs.com/package/dayjs

[annotated functions]: ../../overview/functions/function.md

[applications]: ../how-to/build-an-app.md#applications

[competitor for jQuery]: https://www.npmjs.com/package/cash-dom

[connections]: ../../access/data-connection.md "A data connection is a configuration in Datagrok that lets you access data in a data source such as GitHub repository or local file system."

[custom viewer built with D3.js]: https://github.com/datagrok-ai/public/tree/master/packages/Viewers

[database queries]: ../../access/data-query.md "A data query extract data from a source. A data query can be an SQL query or a query to an Excel file."

[Datagrok configuration section]: _datagrok-configuration.md

[Datagrok tools section]: https://github.com/datagrok-ai/public/tree/master/tools#datagrok-tools

[debug]: ./_debugging.md

[deployment with Docker Compose section]: https://datagrok.ai/help/develop/admin/docker-compose

[development servers]: https://dev.datagrok.ai/u

[functions]: ../../overview/functions/function.md "A function can be used to query a database, delete a column from a table, applying predictive model to a dataset and other things."

[function types]: ./_package-function-types.md

[grok CLI section]: https://github.com/datagrok-ai/public/tree/master/tools#datagrok-tools

[header parameters]: ../../compute/scripting.md#header-parameters

[javaScript API]: https://datagrok.ai/help/develop/js-api

[LTS version of Node.js]: https://nodejs.org/en/download/

[node.js > os homedir section]: https://nodejs.org/api/os.html#os_os_homedir

[npm instructions]: https://docs.npmjs.com/downloading-and-installing-node-js-and-npm

[npm package naming rules]: https://docs.npmjs.com/cli/v6/configuring-npm/package-json#name

[npm]: https://docs.npmjs.com/downloading-and-installing-node-js-and-npm

[nvm]: https://github.com/nvm-sh/nvm#installing-and-updating

[package section]: https://datagrok.ai/js-api/classes/dg.package

[public repository]: https://github.com/datagrok-ai/public/tree/master/packages

[public]: https://public.datagrok.ai/u

[publish the package]: _publishing.md

[script for Cheminformatics]: https://github.com/datagrok-ai/public/tree/master/packages/Chem/scripts

[scripts]: ../../compute/scripting.md "Scripting combines fast interactive visualizations and other features of the Datagrok platform with statistical packages and visualizations available in R, Python, Octave, Julia, and JavaScript."

[semantic type detectors]: ../how-to/define-semantic-type-detectors.md

[sequence package]: https://github.com/datagrok-ai/public/tree/master/packages/Sequence

[tables]: ../../access/connectors/files.md#supported-tabular-formats

[user profile on Datagrok]: https://dev.datagrok.ai/u

[viewers]: ../../visualize/viewers.md "A viewer is a visual component associated with a table."

[views]: ../how-to/custom-views.md

[webpack documentation]: https://webpack.js.org/configuration/

[widgets]: ../../visualize/widgets.md "Widgets are various UI elements that together comprise the platform's user interface."