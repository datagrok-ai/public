<!-- TITLE: JavaScript Development -->
<!-- SUBTITLE: -->

# JavaScript Development

Datagrok was designed to be as extensible as possible, so naturally JavaScript-based development
is the preferred way to develop user-facing applications on top of the platform. 
Use [Grok API](js-api.md) to control pretty much anything within Datagrok,
including [data manipulation](js-api.md#data-manipulation), 
adding [views](js-api.md#views) or [viewers](js-api.md#pre-defined-viewers), 
[developing custom viewers](js-api.md#custom-viewers),
[registering functions](js-api.md#registering-functions),
training and applying [predictive models](../learn/predictive-modeling.md), 
and even [building custom apps](app.md).

There are two options to run custom JavaScript code. For ad-hoc scripts, use the built-in
JavaScript editor (`Functions | Scripts | New JavaScript Script`). For reusable functions, viewers, 
and applications, use the packaging mechanism, which is the focus of this article.

Table of contents

* [Packages](#packages)
* [Getting Started](#getting-started)
* [Package Structure](#package-structure)
* [Development](#development)
* [Publishing](#publishing)
* [Applications](#applications)
* [Documentation](#documentation)

## Packages

A package is a versionable unit of content distribution within Datagrok. Essentially, it is 
a folder with files in it. A package might contain different things:

* JavaScript functions, viewers, widgets, applications
* Scripts written in R, Python, Octave, Grok, Julia, JavaScript, NodeJS, or Java
* Queries and connections
* Tables

See our [GitHub repository](https://github.com/datagrok-ai/public/tree/master/packages) for examples.

## Getting Started

To develop a package on the Datagrok platform, you will need [Node.js](https://nodejs.org/en/) and [npm](https://www.npmjs.com/get-npm) installed. Also, install [Webpack](https://webpack.js.org/guides/installation/) to be able to build your package locally and debug it using `Webpack DevServer`. Optionally, you can use `Babel`, `React` as well as other advanced JavaScript frameworks.

Here are the first steps to get you started:

1. Install `datagrok-tools` utility for managing packages:
   ```
   npm install datagrok-tools -g
   ```
2. Configure your environment with the following command:
   ```
   grok config
   ```
   Enter developer keys and set the default server. Your credentials will be stored locally in `config.yaml`. Once created, this file will be used for publishing all your packages. The developer key can be retrieved by opening your user profile and clicking on `Developer key`. Administrators can manage existing keys and grant or revoke privileges.
3. Create a new package by running this command:
   ```
   grok create <packageName>
   ```
   A new folder `MyPackage` will be created automatically as well as its contents.
4. Once you've completed the work on your package, upload it by running:
   ```
   grok publish
   ```

If you are developing a package using the old template, please run `grok migrate`. This command will convert your scripts in `package.json` and copy keys from `upload.keys.json` to `config.yaml`. Run `grok` for instructions and `grok <command> --help` to get help on a particular command.

## Package Structure

A simplest JavaScript package consists of the following files:

| file                          | description     |
|-------------------------------|-----------------|
| [package.json](#package.json) | metadata        |
| [package.js](#package.js)     | entry point     |
| [detectors.js](#detectors.js) | detectors file  |
| README.md                     | package summary |
| package.png                   | package icon    |

### <a name="package.json"></a>package.json

`package.json` contains metadata, such as name, version, and dependencies: 

```json
{
  "name": "sequence",
  "fullName": "Sequence",
  "version": "0.0.1",
  "description": "Support for DNA sequences",
  "sources": ["ntseq/ntseq.js", "feature-viewer/feature-viewer.bundle.js", "msa/msa.min.gz.js"],
  "dependencies": {
    "datagrok-api": "latest"
  },
  "scripts": {
    "debug-sequence": "grok publish --rebuild",
    "release-sequence": "grok publish --rebuild --release",
    "build-sequence": "webpack",
    "build": "webpack"
  }
}
```

`sources` contains a list of JavaScript or CSS files defined relative to the package's root. Each of these files will be loaded before any function from that package is executed.

The package template first includes only one dependency — `datagrok-api`. You can add more packages to the dependencies list and install them via `npm install`.

The file `package.json` also contains `scripts` for [debugging and publishing your package](#publishing).

### <a name="package.js"></a>package.js

Next, let's take a look at the `src/package.js` file:

```js
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: test
//input: string s
export function test(s) {
    grok.shell.info(_package.webRoot);
}
```

Note that `Datagrok API` modules are already imported. They are also set as external modules, so that `Webpack` will not include them to the output. You can include other libraries or packages, as all of them will be built in the single bundle file.

During the [publishing step](#publishing), the contents of `package.js` get parsed, and functions with the properly formatted
[headers](../compute/scripting.md#header) are registered as Grok [functions](../overview/functions/function.md). By annotating
functions in a specific way, it is possible to register custom viewers, widgets, renderers, 
converters, validators, suggestions, info panels, and semantic type detectors. 

### <a name="detectors.js"></a>detectors.js

`detectors.js` is a JavaScript file.  It should define a class named `<package_name>PackageDetectors` that subclasses `DG.Package`.
It is similar to `package.js` but intended for smaller functions — semantic type detectors. Datagrok calls these functions each time the
user opens a table. Detectors will be uploaded separately from the rest of the package and used to quickly inspect the data and determine the semantic type of the columns. Semantic type tagging allows the platform to offer specific functions for data of a particular type.

Below, there is an example of a package `Sequence` containing a single detector `detectNucleotides`:

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

Once registered, this function is now available across the whole platform, and can be used for semantic type detection.

## Development

A JavaScript package runs inside the Datagrok platform, which makes the development and
debugging experience different compared to the more traditional web applications. Essentially, the packages
are developed locally, but the platform runs remotely. To enable the best possible experience for
developers, we established a workflow where the package is uploaded to the remote server at startup, and
then gets served from the server. By associating local JavaScript files with the remote sources in your favorite IDE, 
it is possible to hide the complexity of that scenario. For instance, you can set breakpoints, do 
step-by-step execution and generally debug the program in the regular way. Of course, you can always
use the debugger that comes with the browser.

To develop Datagrok packages, we recommend that you start with creating a package template.
Then, set up your IDE in such a way that when starting a project, it would [publish](#publishing) the package,
and then start the platform.

Packages deployed in the development mode are visible only to the authors. This ensures that multiple
people can simultaneously work on the same package.    

### General Notes on Package Development

Our approach to extending the system is providing one canonical, convenient way for 
developers to achieve the task, at the same time exposing enough extension points
to enable deep customization. Same concepts apply to JavaScript development. 
We do not impose any requirements on the UI frameworks or technologies used for 
the JavaScript plugins, although we encourage developers to keep it simple.

To simplify development, Datagrok provides an `Inspector` tool (`Alt + I`) that lets developers peek
under the hood of the platform. Use it for understanding which events get fired and when,
how views and viewers are serialized, what is getting stored locally, what widgets are
currently registered by the system, etc. 

### Environments

In order to isolate packages being debugged from the production instance, we recommend running 
them against the `dev` instance, if possible. To change Datagrok's server, add a new developer key to your local `config.yaml` and edit the `scripts` section in the `package.json` file.

## Publishing

### Version Control

Each package has a version number. All objects inside the package are being deployed according to the package version. 
When a package gets published, a "published package" entity gets created. It is associated with the 
package, and has additional metadata (such as publication date). Typically, only one version of a package
is visible to a user. Administrators can manage published packages and decide which versions
should be used. It is possible to roll back to an older version, or assign a particular version to
a particular group of users.

Importantly, if the version changes, there will be an independent instance of each package asset.  
Multiple versions of a package can be deployed at one moment, and the administrator can switch between them. All users 
will only see objects that belong to the current package version.

There is a special `debug` version that can be deployed for each package. If the developer applies it, it becomes 
active for the current package until the developer deletes it or changes their developer key. In this case, the developer can see objects from their version of package, and changes will not affect other users package representation. 

### Deployment Modes

You can use the following flags to specify who can access your package:

* In `--debug` mode, packages are accessible by the developer only (default).
* In `--release` mode, packages are accessible by everyone who has the privilege.

To publish a package, open the `package.json` file in your IDE. Typically, the `scripts` section would contain several scripts generated for your package based on the contents of `config.yaml`. For development purposes, use the scripts having the `debug` word in their name. For production, use the alternative scripts starting with `deploy` instead.

```json
"scripts": {
  "debug-sequence": "grok publish --rebuild",
  "release-sequence": "grok publish --rebuild --release",
  "build-sequence": "webpack",
  "build": "webpack",
  "debug-sequence-dev": "grok publish dev --rebuild",
  "release-sequence-dev": "grok publish dev --rebuild --release",
  "debug-sequence-local": "grok publish local --rebuild",
  "release-sequence-local": "grok publish local --rebuild --release"
}
```

Alternatively, you can run these scripts with `npm run`. If you have any questions, type `grok` for instructions or `grok publish --help` to get help on this particular command.

In addition, you can pass another server either as URL or server alias from the `config.yaml` file:

```
grok publish dev
grok publish https://dev.datagrok.ai/api --key dev-key
```

Make sure to specify the developer key for a new server.

### Source Control

Packages can be deployed from git as well as other resources, which allows for convenient team collaboration and version management. See the full list of source types in the [Package Browser](https://public.datagrok.ai/packages) (`Manage | Packages | Add new package`). When developing a package with your team, it's a good idea to commit code to the repository first and then publish your package from there. Our [public GitHub repository](https://github.com/datagrok-ai/public/tree/master/packages) is a telling example of this workflow. We also welcome contributions, which you can learn more about in [this article](https://datagrok.ai/help/develop/public-repository).

### Continuous Integration

`Webpack` is required for your package source code to work successfully in the browser. The Datagrok platform can build a package on the server side, in which case, you need to specify `--rebuild` option in scripts from `package.json`. The script `"build": "webpack"` is reserved for server build.

Alternatively, it's possible to build your package on the client side using `Webpack`. To do that, specify the opposite option `--build` in your scripts.

### Sharing

Just like other entities on the platform, packages are subject to [privileges](../govern/security.md#privileges). When sharing with users and groups of users, you can specify the rights (for viewing and editing) and choose if you want to notify the person in question.

To see packages available to you, click on `Manage | Packages`, or follow [this link](https://public.datagrok.ai/packages) from outside the platform.

## Applications

Applications are [functions](../overview/functions/function.md) tagged with the `#app` tag. A package might contain zero, one, or more apps. See our [GitHub repository](https://github.com/datagrok-ai/public/tree/master/packages) for application examples, such as [Enamine Store application](https://github.com/datagrok-ai/public/tree/master/packages/EnamineStore).

To open the application launcher, click on `Functions | Apps`, or follow [this link](https://public.datagrok.ai/apps)
from outside the platform. To launch a particular app automatically, open the following URL: `https://public.datagrok.ai/apps/<APP_NAME>`

## Documentation

According to [this study](http://sigdoc.acm.org/wp-content/uploads/2019/01/CDQ18002_Meng_Steinhardt_Schubert.pdf),
in terms of the strategies used for understanding API documentation, different developers fall into several groups:
systematic, opportunistic and pragmatic. These findings are consistent with our experience. 
For Datagrok's documentation, we have established an approach that enables developers from either of the
above-mentioned groups to be productive. 

* [Sample browser](https://public.datagrok.ai/js) (`Functions | Scripts | New JavaScript Script`) is an interactive tool
  for browsing, editing, and running JavaScript samples that come with the platform. Samples are grouped by
  domain, such as data manipulation, visualization, or cheminformatics. They are short, clean examples
  of the working code using [Grok API](js-api.md) that can be copy-and-pasted into the existing solution.
  The samples are also cross-linked with the [help](https://datagrok.ai/help) system. 
* [Grok API](js-api.md) provides complete control over the platform. 
  [JS documentation](https://public.datagrok.ai/js) is available.
* [Platform help](https://datagrok.ai/help/) explains the functionality from the user's point of view. Where
  appropriate, it is hyper-linked to samples and demo projects. In the near future, we plan to turn it
  into the community wiki, where users will be contributing to the content. The same web pages
  are used as an interactive help within the platform (you see help on the currently selected object).
  
Additionally, there are a few ways to connect with fellow developers:    
* [Datagrok community](https://community.datagrok.ai/)
* [Slack space](https://datagrok.slack.com) 

See also: 

* [Grok API](js-api.md)
* [Packages from our GitHub repository](https://github.com/datagrok-ai/public/tree/master/packages)
* [How Developers Use API Documentation: An Observation Study](http://sigdoc.acm.org/wp-content/uploads/2019/01/CDQ18002_Meng_Steinhardt_Schubert.pdf)
