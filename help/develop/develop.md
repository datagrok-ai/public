<!-- TITLE: &#8204;JavaScript development-->
<!-- SUBTITLE: -->

# JavaScript development

JavaScript or TypeScript-based development is the preferred way to develop user-facing applications on top of the
platform. Use the [JS API](js-api.md) to control pretty much anything within Datagrok,
including [data manipulation](js-api.md#data-manipulation), adding [views](js-api.md#views)
or [viewers](how-to/manipulate-viewers.md),
[developing custom viewers](how-to/develop-custom-viewer.md),
[registering functions](js-api.md#registering-functions), training and
applying [predictive models](../learn/predictive-modeling.md), and even [building custom apps](#applications).

There are two options to run custom JavaScript code. For ad-hoc [scripts](../compute/scripting.md), use the built-in
JavaScript editor (`Functions | Scripts | New JavaScript Script`). For reusable functions, viewers, and applications,
use the packaging mechanism, which is the focus of this article.

This article describes what a [package](#packages) is, as well as techniques for [developing](#development),
[debugging](#debugging), [publishing](#publishing) and using [documentation](#documentation).

## Packages

A package is a versionable unit of content distribution within Datagrok. Essentially, it is a folder with files in it. A
package might contain different things:

* JavaScript [functions](../overview/functions/function.md), [viewers](../visualize/viewers.md)
  , [widgets](../visualize/widgets.md), [applications](#applications)
* [Scripts](../compute/scripting.md) written in R, Python, Octave, Grok, Julia, JavaScript, NodeJS, or Java
* [Queries](../access/data-query.md) and [connections](../access/data-connection.md)
* [Tables](../access/connectors/files.md#supported-tabular-formats), files, and other objects

See our [GitHub repository](https://github.com/datagrok-ai/public/tree/master/packages) for examples, or follow
the [step-by-step guide](how-to/create-package.md) for creating your own package.

## Package structure

A simplest JavaScript package consists of the following files:

| file                                    | description           |
|-----------------------------------------|-----------------------|
| [package.json](#package.json)           | metadata              |
| [package.js](#package.js)               | entry point           |
| [detectors.js](#detectors.js)           | detectors file        |
| [webpack.config.js](#webpack.config.js) | webpack configuration |
| README.md                               | package summary       |
| package.PNG                             | package icon          |

In addition to that, it might contain the following folders:

* `environments`: [environment configurations](../compute/scripting.md#environments)
  for [scripts](../compute/scripting.md). Examples: [DemoScripts]
* `scripts`: a collection of [scripts](../compute/scripting.md) used for computations.
  Examples: [Chem](https://github.com/datagrok-ai/public/tree/master/packages/Chem)
  , [DemoScripts]
  , [Impute](https://github.com/datagrok-ai/public/tree/master/packages/Impute)
* `swaggers`: REST APIs in [Swagger/OpenAPI](../access/open-api.md) format.
  Examples: [EnamineStore](https://github.com/datagrok-ai/public/tree/master/packages/EnamineStore)
  , [Swaggers](https://github.com/datagrok-ai/public/tree/master/packages/Swaggers)
* `connections` and `queries`: [connections](../access/data-connection.md)
  and [queries](../access/data-query.md) for data retrieval.
  Examples: [Chembl](https://github.com/datagrok-ai/public/tree/master/packages/Chembl)
  , [UsageAnalysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
* `css`: CSS files for custom styling.
  Examples: [Notebooks](https://github.com/datagrok-ai/public/tree/master/packages/Notebooks)
  , [Discovery](https://github.com/datagrok-ai/public/tree/master/packages/Discovery)
* `tables` and `data-samples`: data for demonstration and testing.
  Examples: [Chem](https://github.com/datagrok-ai/public/tree/master/packages/Chem)
  , [Sunburst](https://github.com/datagrok-ai/public/tree/master/packages/Sunburst)
* `layouts`: `json` files with table view [layouts](how-to/layouts.md)
* `schemas`: `yaml` files with property schemas
* `jobs`: [data jobs](../access/data-job.md)

### <a href="#" id="package.json"></a>package.json

`package.json` contains metadata, such as name, version, and dependencies:

```json
{
  "name": "sequence",
  "fullName": "Sequence",
  "version": "0.0.1",
  "description": "Support for DNA sequences",
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

The package template first includes only one dependency — `datagrok-api`. You can add more packages to the dependencies
list and install them via `npm install`.

The file `package.json` also contains `scripts`
for [debugging and publishing your package](#publishing).

### <a href="#" id="package.js"></a>package.js

Next, let's take a look at the `src/package.js` file:

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

Note that `Datagrok API` modules are already imported. They are also set as external modules, so that `Webpack` will not
include them to the output. You can include other libraries or packages, as all of them will be built-in a single bundle
file. If you choose to include other files, such as CSS, in your package, import them into `package.js` as well.

During the [publishing step](#publishing), the contents of `package.js` get parsed, and functions with the properly
formatted
[headers](../compute/scripting.md#header) are registered as Grok [functions](../overview/functions/function.md)
. By annotating functions in a specific way, it is possible to register custom viewers, widgets, renderers, converters,
validators, suggestions, info panels, and semantic type detectors. If function has more than one output, it must return
JS object `{param1: value, param2: value}`:

```js
//name: test
//output: string s1
//output: string s1
export function test() {
  return {s1: 'a', s2: 'b'};
}
```

### <a href="#" id="detectors.js"></a>detectors.js

`detectors.js` is a JavaScript file. It should define a class named `<package_name>PackageDetectors`
that subclasses `DG.Package`. It is similar to `package.js` but intended for smaller functions — semantic type
detectors. Datagrok calls these functions each time the user opens a table. Detectors will be uploaded separately from
the rest of the package and used to quickly inspect the data and determine the semantic type of the columns. Semantic
type tagging allows the platform to offer specific functions for data of a particular type.

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

### <a href="#" id="webpack.config.js"></a>webpack.config.js

The package is built according to its configuration file, which typically has the following content:

```javascript
const path = require("path");

module.exports = {
  mode: "development",  // set to "production" to minify the output and enable optimizations for production builds
  entry: {
    package: "./src/package.js"  // the package is limited to exactly one entry point
  },
  devtool: "inline-source-map",   // enhances package debugging in the browser devtools
  externals: {                    // external modules won't be loaded to the output, but taken from the environment
    "datagrok-api/dg": "DG",
    "datagrok-api/grok": "grok",
    "datagrok-api/ui": "ui",
    "openchemlib/full.js": "OCL",
    "rxjs": "rxjs",
    "rxjs/operators": "rxjs.operators"
  },
  output: {
    filename: "[name].js",
    library: "sequence",     // the name of the package in lower case
    libraryTarget: "var",    // the results will be assigned to a variable `sequence`
    path: path.resolve(__dirname, "dist"),
  },
};
```

Have a look at the [Webpack documentation](https://webpack.js.org/configuration/) in case you need to modify or extend
the provided options. For instance, you can add CSS and other file [loaders](https://webpack.js.org/loaders/)
to `module.rules`. When the package is loaded, the output gets assigned to a variable (
type `window.<package_name>`, e.g. `window.sequence`, in the browser's console just to check). Finally, note that the
package name have reoccurred in multiple files, including this one. This might become important if you are going to
introduce changes to the code or, for example, rename the package without creating it from scratch. In this case, make
sure the name is accurately substituted: set the `name` field in `package.json` and `library` in `webpack.config.js` to
the desired name in lower case, and rename a class `<package_name>PackageDetectors` using camel case in `detectors.js`.

### Naming conventions

Continuing the topic we have just touched on, here are naming guidelines and general recommendations that you might
consider:

* Use upper camel case for package names, for example, `ApiSamples` and `OctaveScripts`. Package names that comply with
  the [rules](https://docs.npmjs.com/cli/v6/configuring-npm/package-json#name) for `npm` packages, e.g. `api-samples`
  and `octave-scripts`, are accepted as well. That being said, you can still write the desired name in the `fullName`
  field of `package.json`
* When defining new [views](how-to/custom-views.md) and [viewers](how-to/develop-custom-viewer.md), we recommend
  postfixing your classes with `View` and `Viewer` respectively
* The names of semantic type detectors typically start with the `detect` prefix, e.g., `detectNucleotides`
  or `detectRDSmiles`
* Filenames can be written in lower case, with dashes between words: `tika-extractor.py`
  and `chord-viewer.js`

## Development

A JavaScript package runs inside the Datagrok platform, which makes the development and debugging experience different
compared to the more traditional web applications. Essentially, the packages are developed locally, but the platform
runs remotely. To enable the best possible experience for developers, we established a workflow where the package is
uploaded to the remote server at startup, and then gets served from the server. By associating local JavaScript files
with the remote sources in your favorite IDE, it is possible to hide the complexity of that scenario. For instance, you
can set breakpoints, do step-by-step execution and generally debug the program in the regular way. Of course, you can
always use the debugger that comes with the browser.

To develop Datagrok packages, we recommend that you start with creating a package template. Then, set up your IDE in
such a way that when starting a project, it would [publish](#publishing) the package, and then start the platform.

Packages deployed in the development mode are visible only to the authors. This ensures that multiple people can
simultaneously work on the same package.

### General notes on package development

Our approach to extending the system is providing one canonical, convenient way for developers to achieve the task, at
the same time exposing enough extension points to enable deep customization. Same concepts apply to JavaScript
development. We do not impose any requirements on the UI frameworks or technologies used for the JavaScript plugins,
although we encourage developers to keep it simple.

To simplify development, Datagrok provides an `Inspector` tool (`Alt + I`) that lets developers peek under the hood of
the platform. Use it for understanding which events get fired and when, how views and viewers are serialized, what is
getting stored locally, what widgets are currently registered by the system, etc.

### Environments

In order to isolate packages being debugged from the production instance, we recommend running them against the `dev`
instance, if possible. To change Datagrok's server, add a new developer key to your local `config.yaml` and edit
the `scripts` section in the `package.json` file.

## Publishing

### Version control

Each package has a version number. All objects inside the package are being deployed according to the package version.
When a package gets published, a "published package" entity gets created. It is associated with the package, and has
additional metadata (such as publication date). Typically, only one version of a package is visible to a user.
Administrators can manage published packages and decide which versions should be used. It is possible to roll back to an
older version, or assign a particular version to a particular group of users.

Importantly, if the version changes, there will be an independent instance of each package asset. Multiple versions of a
package can be deployed at one moment, and the administrator can switch between them. All users will only see objects
that belong to the current package version.

There is a special `debug` version that can be deployed for each package. If the developer applies it, it becomes active
for the current package until the developer deletes it or changes their developer key. In this case, the developer can
see objects from their version of package, and changes will not affect other users package representation. This version
will no longer exist after the developer releases their package.

### Deployment modes

Use the following flags to specify who can access your package:

* In `--debug` mode, packages are accessible by the developer only (default).
* In `--release` mode, packages are accessible by everyone who has the privilege.

To publish a package, open the `package.json` file in your IDE. Typically, the `scripts` section would contain several
scripts generated for your package based on the contents of `config.yaml`. For development purposes, use the scripts
having the `debug` word in their name. For production, use the alternative scripts starting with `deploy` instead.

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

Alternatively, you can run these scripts with `npm run`. If you have any questions, type `grok` for instructions
or `grok publish --help` to get help on this particular command.

In addition, you can pass another server either as URL or server alias from the `config.yaml` file:

```
grok publish dev
grok publish https://dev.datagrok.ai/api --key <dev-key>
```

Make sure to specify the developer key for a new server.

### Source control

Packages can be deployed from `Git` as well as other resources, which allows for convenient team collaboration and
version management. See the full list of source types in
the [Package Browser](https://public.datagrok.ai/packages) (`Manage | Packages | Add new package`).

When developing a package with your team, it's a good idea to commit code to the repository first and then publish your
package from there. Our [public GitHub repository](https://github.com/datagrok-ai/public/tree/master/packages) is a
telling example of this workflow. We also welcome contributions, which you can learn more about
in [this article](https://datagrok.ai/help/collaborate/public-repository).

To publish a package from the repository, you need to open `Manage | Packages | Add new package`
first. Once the window appears, choose `Git` as the source type, enter the URL to your repository, and specify the
package directory relative to its root. Click on `LOAD PACKAGE METADATA` to get the package name and description.

![Git Publish](img/git-publishing.png)

If necessary, you can specify additional settings and then publish the package.

### Continuous integration

`Webpack` is required for your package source code to work successfully in the browser. The Datagrok platform can build
a package on the server side, in which case, you need to specify `--rebuild`
option in scripts from `package.json`. The script `"build": "webpack"` is reserved for server build.

Alternatively, it's possible to build your package on the client-side using `Webpack`. To do that, run the
script `"build": "webpack"` before publishing.

Package publication is compatible with automation tools. You can pass your server URL and developer key explicitly
without configuring:

```
grok publish <url> -k <dev-key>
```

### Sharing

Just like other entities on the platform, packages are subject to [privileges](../govern/security.md#privileges). When
sharing with users and groups of users, you can specify the rights (for viewing and editing) and choose if you want to
notify the person in question. These privileges can be managed not only from the user interface, but also directly from
the package. To do that, you should specify the eligible user groups in the `package.json` file:

```json
{
  "canEdit": [
    "Developers"
  ],
  "canView": [
    "All users"
  ]
}
```

To see packages available to you, click on `Manage | Packages`, or
follow [this link](https://public.datagrok.ai/packages) from outside the platform.

## Debugging

Essentially, Datagrok packages are either regular webpack packages, or sets JS source files.

We recommend you the practices below with debugging your package code.

### One-click debugging with Visual Studio Code

Use an `--ide=vscode` key of [grok tools](#Getting-Started) to set up a `.vscode` debugging integration with Datagrok.

This parameter coupled with `grok create` adds an additional `.vscode` folder with two files to your package folder.
These files bring you development experience close to building native apps, where you modify the code, hit the "Run",
wait for application to "compile" (in our case — go through webpack, and potentially through linters), start it and then
explore your runtime state through breakpoints, variable watches, call stacks etc.

With this setting, you would be able to "compile" and "run" the package in one click from the left sidebar of Visual
Studio Code under the `Run` section (4-th button from the top) by the `Start Debugging` button at the top. The steps
added with grok tools `--vscode` key will prepare a webpack package locally (showing occurring errors accordingly),
deploy it to the selected Datagrok instance, and run the browser in a IDE-debugging mode.

Below we explain how this debugging works under the hood. As of today, you'd also need to repeat these steps for VS Code
on Linux and MacOS, as we don't support automatically created configurations for them yet.

The first file, `launch.json`, has contents similar to the below:

```
{
  "version": "0.2.0",
  "configurations": [
    {
      "preLaunchTask": "rebuild",
      "type": "chrome",
      "request": "launch",
      "name": "Launch Chrome against Dev",
      "url": "https://dev.datagrok.ai/",
      "webRoot": "${workspaceFolder}"
    }
  ]
}
```

It allows VS Code to launch an isolated Chrome session in debugging mode linked to IDE, targeting the desired
web-address from `url`. In case an application is being developed, it may be useful to target the application's URL.

The second file, `tasks.json`, contains the following (for Windows):

```
{
  "version": "2.0.0",
  "tasks": [
    {
      "type": "shell",
      "command": "cmd.exe /c \"webpack && grok publish dev\"",
      "label": "rebuild"
    }
  ]
}
```

This command builds your package with webpack. The webpack step helps track some syntactic errors before they occur in
runtime, and consequently ship the package to your Datagrok server. In this example, a default server would be chosen to
deploy to. The server configurations are pre-created by you with `grok config` eariler and stored
in `%USERPROFILE%/.grok/config.yaml`, where you can modify them, and choose the default configuration. To override the
default configuration, call `grok publish`
with an explicit configuration name: i.e. `grok publish dev`.

After these two files are provided in the working folder, what is left is to set up VS Code internally.

* Open `Run → Install additional debuggers...`, find `Debugger for Chrome` in the appeared list and install it; restart
  VS Code.

* Activate `View → Appearance → Show activity bar`. This will bring the `Run` button into the view. You can also
  use `Ctrl+Shift+D` to run the application.

Make sure the `.vscode` folder is inside the folder you open with VS Code via `File → Open Folder`.

Now open the folder with your package with `File → Open Folder`, add some breakpoint to your code, and hit the `Run`
button. You should see a view with a left-side pane and a green "Run" (play) button on top of it, with the neightbour
dropdown repeating the name of the configuration from `launch.json`. If this is not the case, and you see some other
names or an invitation to create a configuration, it means the pre-made `.vscode` folder isn't recognized by VS Code. In
such case, check that you do `File → Open Folder` on the folder which *contains* the pre-made `.vscode` folder.

First time you hit run, you'll need to enter your Datagrok credentials into the newly appeared Chrome window. After you
have done this, close Chrome, stop debugging from VS Code and start debugging same way again. This time your breakpoints
should work, and you'd see locals and stack traces as usual.

This [video](https://youtu.be/zVVmlRorpjg?list=PLIRnAn2pMh3kvsE5apYXqX0I9bk257_eY&t=871) gives an overview of setting up
the debugging in VS Code.

### One-click debugging with JetBrains IDEs

Below we use WebStorm IDE on Windows as an example, this practice should extrapolate to IntelliJ IDEA, and to other OS.

Open the `Edit Configurations...` dialog in WebStorm:

![WebStorm: Edit Configurations](packages/webstorm-debugging-01.png)

In this dialog, add a new `Shell Script` configuration by the `+` button:

![WebStorm: Adding a Shell Script configuration](packages/webstorm-debugging-02.png)

Fill this configuration with the following:

![WebStorm: Shell Script configuration content](packages/webstorm-debugging-03.png)

*Note.* This configuration looks controversial, but this is so only to overcome the known problems in WebStorm. It turns
out both interpreter and script parameters have to be specified. Omitting the script parameter will work if you run
the `Webpack & Publish` configuration standalone, but *won't*
work as part of the other's configuration build step, which we'd need later. As for the trailing `echo`, it is simply
due to the way this build step works in WebStorm:
it concatenates the interpreter line with the script call line, and we use it to suppress the redundant script call
line.

This `Webpack & Publish` configuration shall become a step running a webpack on your package and publishing it to a
configuration of your choice. In the above screenshot we've chosen a `dev`
configuration, declared in your grok config file at `%USERPROFILE%/.grok/config.yaml`:

![grok config settings](packages/webstorm-debugging-04.png)

Use the configuration of your choice in the `Shell Script` configuration instead of `dev`.

Proceed to creating the actual JavaScript debug configuration. Using the same `+` button, add a `JavaScript Debug`
configuration, and fill it with the following:

![WebStorm: JavaScript Debug configuration content](packages/webstorm-debugging-05.png)

Basically, you need to add the previously created `Shell Script` configuration as a `Before Launch`
step. After this step is executed on hitting the `Debug` button and displaying the `webpack` log and errors in the
WebStorm tool output window, WebStorm will launch the browser in the debugging mode and let you hit your breakpoints:

![WebStorm: a debugging session](packages/webstorm-debugging-06.png)

*Troubleshooting.* Many users reported inability to employ any kind of JavaScript debugging from the WebStorm IDE after
it initially worked. As of 2021, this is
a [known issue](https://intellij-support.jetbrains.com/hc/en-us/community/posts/360009567459-Webstorm-2020-2-1-Remote-Debugging-do-not-work)
of WebStorm IDEs. To fix it, remove these two files and restart the IDE (replace the IDE version with yours):

* `%USERPROFILE%\AppData\Roaming\JetBrains\WebStorm2020.1\options\web-browsers.xml`
* `%USERPROFILE%\AppData\Roaming\JetBrains\WebStorm2020.1\options\other.xml`

*Note 1.* You may notice there's an option to add a `Run External tool` step, which could serve us
this `Webpack & Publish`
step. Unfortunately, this step won't work due to a [known issue](https://youtrack.jetbrains.com/issue/IDEA-229467)
in WebStorm IDEs. We broke the run step into two steps intentionally to alleviate this problem.

*Note 2.* Since 2017 the JetBrains IDE Support plugin is
[no longer required](https://intellij-support.jetbrains.com/hc/en-us/community/posts/360010507240-where-is-JETBRAINS-IDE-SUPPORT-chrome-extension-it-cant-be-found-anywhere-now-on-the-internet)
to debug JavaScript from JetBrains IDEs. In February 2021 JetBrains has removed this plugin from the Chrome store.

*See also:*

* Debugging JavaScript in WebStorm ([1](https://www.youtube.com/watch?v=Qcqnmle6Wu8),
  [2](https://www.youtube.com/watch?v=YNNDMpoGV0w))
* IntelliJ IJEA JavaScript debugging ([link](https://www.jetbrains.com/help/idea/debugging-javascript-in-chrome.html))

### Debugging through the browser

Open Chrome Developers Console (hit `F12`) and go to `Sources` tab.

#### Webpack-based packages

If you deploy package in debug mode (`--release` isn't passed to `grok publish`) and source maps are properly
generated (thus, the setting [devtool](https://webpack.js.org/configuration/devtool/) of
`module.exports =` section in `webpack.config.js` is present), you'd find your package sources in the `top` (root)
section of the source tree by its decapitalized name.

#### Source-based packages

Deploying such package locates it to the Datagrok host URI (such as `https://dev.datagrok.ai`) under
`api → packages/published/flies → <PACKAGE_NAME>/<VERSION>/_/<DIGIT>`, where you'd set breakpoints.

### Troubleshooting debugging

1. For webpack-based packages, make sure there is `devtool: 'inline-source-map'`
   in `module.exports =` section of
   `webpack.config.js`. Otherwise, source maps aren't generated and the IDE won't get source code locations.

2. Make sure the required plugins / debuggers for Chrome debugging are installed in your IDE.

## Package settings

A package can have settings, which are either set programmatically or by users in the package's properties panel. Every
user group has its own settings configuration. In the interface, users will be able to adjust the settings for each
group they belong to. To include a settings editor into a package, add the list of properties to the `package.json`
file:

```json
"properties": [
{
"name": "Property name",
"propertyType": "string", // `DG.TYPES_SCALAR` are supported
"choices": ["value #1", "value #2"], // Optional field with values of the property type
"defaultValue": "value #2", // Optional field with a default value (it should in choices, if they are given)
"nullable": false                     // Optional field determining whether the property value can be null
}
]
```

To retrieve the state of package settings in code, use the `getProperties` method of `DG.Package`:

```js
const props = await _package.getProperties();
```

The above call outputs an object where the keys are property names and the values are serialized property values. It's
possible to customize the editor's appearance by defining a special [editor function](#settings-editors).

## Documentation

According to [this study](http://sigdoc.acm.org/wp-content/uploads/2019/01/CDQ18002_Meng_Steinhardt_Schubert.pdf)
, in terms of the strategies used for understanding API documentation, different developers fall into several groups:
systematic, opportunistic and pragmatic. These findings are consistent with our experience. For Datagrok's
documentation, we have established an approach that enables developers from either of the above-mentioned groups to be
productive.

* [Sample browser](https://public.datagrok.ai/js) (`Functions | Scripts | New JavaScript Script`) is an interactive tool
  for browsing, editing, and running JavaScript samples that come with the platform. Samples are grouped by domain, such
  as data manipulation, visualization, or cheminformatics. They are short, clean examples of the working code
  using [Grok API](js-api.md)
  that can be copy-and-pasted into the existing solution. The samples are also cross-linked with
  the [help](https://datagrok.ai/help) system.
* [Grok API](js-api.md) provides complete control over the platform.
  [JS documentation](https://public.datagrok.ai/js) is available.
* [Platform help](https://datagrok.ai/help/) explains the functionality from the user's point of view. Where
  appropriate, it is hyper-linked to samples and demo projects. In the near future, we plan to turn it into the
  community wiki, where users will be contributing to the content. The same web pages are used as an interactive help
  within the platform (you see help on the currently selected object).

Also, you can connect with fellow developers on either
[community forum](https://community.datagrok.ai/) or [slack](https://datagrok.slack.com).

See also:

* [Grok API](js-api.md)
* [Scripting](../compute/scripting.md)
* [Packages from our GitHub repository](https://github.com/datagrok-ai/public/tree/master/packages)
* [How Developers Use API Documentation: An Observation Study](http://sigdoc.acm.org/wp-content/uploads/2019/01/CDQ18002_Meng_Steinhardt_Schubert.pdf)

[DemoScripts]: https://github.com/datagrok-ai/public/blob/master/packages/Demo/projects/scripts/
