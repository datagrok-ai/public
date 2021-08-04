# Packages

You can extend Datagrok with _packages_ that you build using the [JavaScript API]. A Datagrok package is a JavaScript or 
TypeScript application that introduces new functionality to Datagrok or extends the existing functionality.

For example, you package might create [custom viewers using D3], demonstrate [scripting functionality for Cheminformatics],
and do even more. Check out the [public repository] for more examples of Datagrok packages.

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
* [Tables]

## Creating a package

To develop a package on the Datagrok platform, install the [latest LTS version of Node.js] and [npm]. We strongly
recommend that you install Node.js using [NVM] to avoid permission issues when running NPM packages globally. Check out
the [NPM instructions] for more details.

Follow these steps to create a package template:

1. Install `datagrok-tools` globally to manage your Datagrok packages:

```shell
npm install datagrok-tools -g
```

2. Configure your environment:

```shell
grok config
```

You will be prompted to enter the developer keys and set the default server. Your credentials will be stored locally
in `config.yaml`. Check out the [Datagrok configuration section] for more details.

3. Go to the folder where you want to create your package and run:

```shell
grok create <package-name>
```

> **Note**: If you want to create a TypeScript package, pass the `--ts` option to the command.

After you complete the work on your Datagrok package, upload it to Datagrok by running the following command:

```shell
grok publish
```

Note that by running this command, only you will see the package in Datagrok. For more information about the available
commands, refer to the [Grok CLI section].

## Package structure

| File                                    | Description                                               |
|-----------------------------------------|---------------------------------------------------------- |
| node_modules                            | Node.js modules                                           |
| [src/package.js](#package.js)           | The Datagrok package's entry point                        |
| .gitignore                              | Files and folders ignored by Git                          |
| [detectors.js](#detectors)              | File with detectors                                       |
| [package.json](#package.json)           | Package metadata with predefined dependencies and scripts |
| package.png                             | Package icon                                              |
| README.md                               | Summary of the package                                    |
| [webpack.config.js](#webpack.config.js) | webpack configuration for development                     |

We discuss the key files in detail in their own sections that follow.

### <a href="#" id="package.json"></a>package.json

The `package.json` file contains the project metadata.

```json
{
  "name": "test-package",
  "beta": true,
  "friendlyName": "test-package",
  "version": "0.0.1",
  "description": "",
  "dependencies": {
    "datagrok-api": "latest"
  },
  "devDependencies": {
    "webpack": "latest",
    "webpack-cli": "latest"
  },
  "scripts": {
    "debug-test-package": "grok publish --rebuild",
    "release-test-package": "grok publish --rebuild --release",
    "build-test-package": "webpack",
    "build": "webpack",
    "debug-test-package-dev": "grok publish dev --rebuild",
    "release-test-package-dev": "grok publish dev --rebuild --release",
    "debug-test-package-local": "grok publish local --rebuild",
    "release-test-package-local": "grok publish local --rebuild --release"
  }
}
```

The package includes only the `datagrok-api` dependency &mdash; this is the JavaScript API you will actively use when
developing the package. Additionally, the CLI includes `webpack` and `webpack-cli` as development dependencies to let 
you build your project.

You can install other dependencies (such as React) as you need.

`package.json` also contains [a few valuable scripts] to let you [debug](4__debugging.md) or 
[publish the package](5__publishing.md).

### <a href="#" id="package.js"></a>src/package.js

Here's the `src/package.js` file (or `src/package.ts` if you created a TypeScript package):

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

The `package.js` imports three key Datagrok modules &mdash; `datagrok-api/grok`, `datagrok-api/ui`, and 
`datagrok-api/dg`.

This file typically contains the [annotated functions] that implement a functionality. Each function needs to be 
annotated with respective [header parameters] that define the metadata of the function. The parameters help Datagrok to 
register your viewers, widgets, renderers, converters, validators, suggestions, information panels, and semantic type 
detectors.

If a function has more than one output, it must return a JavaScript object `{prop1: value1, prop2: value2}`, for 
example:

```js
//name: test1
//output: string a1
export function test1() {
  return 'a';
}
//name: test2
//output: string s1
//output: string s1
export function test2() {
  return {s1: 'a', s2: 'b'};
}
```

### <a href="#" id="detectors.js"></a>detectors.js

The `detectors.js` file must define a `<packageName>PackageDetectors` class that extends the base `DG.Package` class.
(For more information about this class, refer to the [JavaScript API - Package section].)
`detectors.js` is similar to `package.js` but is intended for smaller functions called [_semantic type detectors_],
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
    mode: 'development',  // Set to "production" to minify the output and optimize production builds
    entry: {
        package: './src/package.js'  // The package is limited to exactly one entry point
    },
    devtool: 'inline-source-map',   // Enhances package debugging in browser devtools
    externals: {                    // The external modules won't be loaded to the output, but taken from the environment
        'datagrok-api/dg': 'DG',
        'datagrok-api/grok': 'grok',
        'datagrok-api/ui': 'ui',
        'openchemlib/full.js': 'OCL',
        'rxjs': 'rxjs',
        'rxjs/operators': 'rxjs.operators'
    },
    output: {
        filename: '[name].js',
        library: 'sequence',     // Name of the package in lower case
        libraryTarget: 'var',    // Results will be assigned to a variable sequence`
        path: path.resolve(__dirname, 'dist'),
    },
};
```

Refer to the [webpack documentation] to modify or extend the provided configuration.

## Naming conventions

Follow these naming conventions when working on your packages:

* Use upper camel case for package names: `ApiSamples` and `OctaveScripts`
  > **Note**: Package names that comply with the [NPM package naming rules], e.g., `api-samples` and `octave-scripts`,
  are accepted as well. If you write the desired name in the `friendlyName` field of `package.json`, it will be shown in the UI.
* When defining new [views] and [viewers], postfix the classes with `View` and `Viewer` respectively
* Prefix the names of semantic type detectors with `detect`: `detectNucleotides` or `detectRDSmiles`
* Use lowercase letters for files, separate words by dashes: `tika-extractor.py` and `chord-viewer.js`

## Structuring your package

Your package might contain the following additional folders, depending on your needs:

| Folder       | Description                            | Documentation               | Package Examples                    |
|------------- | -------------------------------------- | --------------------------- | ----------------------------------- |
| environments | Environment configurations             | [Scripting > Environments]  | [PythonScripts]                     |
| scripts      | Collection of scripts for computations | [Scripting]                 | [ChemScripts], [RScripts], [Impute] |
| swaggers     | OpenAPI specifications                 | [Access > OpenAPI]          | [EnamineStore], [Swaggers]          |
| connections  | Database connections                   | [Access > Data Connections] | [Chembl], [UsageAnalysis]           |
| queries      | Database queries                       | [Access > Data Queries]     | [Chembl], [UsageAnalysis]           |
| css          | CSS files with custom styles           | -                           | [Notebooks], [Discovery]            |
| data-samples | Data for demo and testing              | -                           | [Chem], [Sunburst]                  |

[JavaScript API]: https://datagrok.ai/help/develop/js-api
[custom viewers using D3]: https://github.com/datagrok-ai/public/tree/master/packages/Viewers
[scripting functionality for Cheminformatics]: https://github.com/datagrok-ai/public/tree/master/packages/ChemScripts
[public repository]: https://github.com/datagrok-ai/public/tree/master/packages
[Functions]: ../../overview/functions/function.md
[views]: ../how-to/custom-views.md
[viewers]: ../../visualize/viewers.md
[widgets]: ../../visualize/widgets.md
[applications]: 3__package-function-types.md#applications
[Scripts]: ../scripting.md
[Database queries]: ../../access/data-query.md
[connections]: ../../access/data-connection.md
[Tables]: ../../access/connectors/files.md#supported-tabular-formats
[latest LTS version of Node.js]: https://nodejs.org/en/download/
[npm]: https://docs.npmjs.com/downloading-and-installing-node-js-and-npm
[NVM]: https://github.com/nvm-sh/nvm#installing-and-updating
[NPM instructions]: https://docs.npmjs.com/downloading-and-installing-node-js-and-npm
[NPM package naming rules]: https://docs.npmjs.com/cli/v6/configuring-npm/package-json#name
[Datagrok configuration section]: ./datagrok-configuration.md
[Grok CLI section]: https://github.com/datagrok-ai/public/tree/master/tools#datagrok-tools
[a few valuable scripts]: https://github.com/datagrok-ai/public/tree/master/tools#commands
[annotated functions]: ../../overview/functions/function.md
[header parameters]: ../scripting.md#header-parameters
[JavaScript API - Package section]: https://datagrok.ai/js-api/Package
[_semantic type detectors_]: ../how-to/semantic-type-detector.md
[Sequence package]: https://github.com/datagrok-ai/public/tree/master/packages/Sequence
[webpack documentation]: https://webpack.js.org/configuration/
