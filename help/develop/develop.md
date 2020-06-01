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
JavaScript editor (`Tools | Scripting | JavaScript`). For reusable functions, viewers, 
and applications, use the packaging mechanism.

Table of contents

* [Getting Started](#getting-started)
* [Packages](#packages)
* [Buildable Packages](#buildable-packages)
* [Development](#development)
* [Publishing](#publishing)
* [Documentation](#documentation)
* [Roadmap](#roadmap)

## Getting Started

Here are the typical steps for creating and debugging applications on the Datagrok platform:

1. Download package template: [Open Packages Browser](https://public.datagrok.ai/packages), click on toolbox 
   "Actions | Download package template..." that contains sources, upload-debug.cmd and upload-deploy.cmd
2. Setup development environment (see [video](https://www.youtube.com/watch?v=PDcXLMsu6UM))
3. Debug using [upload-debug.cmd](#publishing) pre-run step
4. Publish using [upload-deploy.cmd](#publishing) 

If you want to build an advanced package, using WepPack builder, or use Babel, React, 
or some other advanced Javascript Frameworks or packages, please, proceed to [Buildable Packages](#buildable-packages) section. 

## Packages

A package is a versionable unit of distribution of content within Datagrok. Essentially, it is 
a folder with files in it. A package might contain different things:

* JavaScript functions, viewers, widgets, applications.
* Scripts written in R, Python, JavaScript, or Java 
* Queries or connections
* Tables

See our [github repository](https://github.com/datagrok-ai/public/tree/master/packages) for examples.

A simplest JavaScript package consists of the following files:

| file          | description    |
|---------------|----------------|
| package.json  | metadata       |
| package.js    | entry point    |
| detectors.js  | detectors file |
| package.png   | icon           |

`package.json` contains metadata, such as name, version, and dependencies: 

```json
{
  "name": "Sequence",
  "version": "0.1",
  "fullName": "Sequence",
  "description": "Support for DNA sequences",
  "sources": ["ntseq/ntseq.js", "feature-viewer/feature-viewer.bundle.js", "msa/msa.min.gz.js"]
}
```

"sources" contains a list of JS or CSS files, defined relative to the package's root, that should be
loaded before any function from that package is executed.

`package.js` is a JavaScript file. It should define a class named <package_name>Package that subclasses GrokPackage.
During the [publishing step](#publishing), its content gets parsed, and functions with the properly formatted
[headers](../compute/scripting.md#header) get registered as Grok [functions](../overview/functions/function.md). By annotating
functions in a specific way, it is possible to register custom viewers, widgets, renderers, 
converters, validators, suggestions, info panels, and semantic type detectors. 

`detectors.js` is a JavaScript file.  It should define a class named <package_name>DetectorsPackage that subclasses GrokPackage.
It's similar to package.js, but intended for small functions - semantic types detectors. Datagrok calls these functions each time
user opens a table.

Below is an example of a package consisting of a single function `complement` that returns a complementary
sequence for the specified nucleotide sequence:

```js
class SequencePackage extends DG.Package {
    
    //description: returns complementary sequence
    //tags: bioinformatics, converter
    //input: string nucleotides {semType: nucleotides}
    //output: string result {semType: nucleotides}
    complement(nucleotides) {
        var seq = new Nt.Seq();
        seq.read(nucleotides);
        return seq.complement().sequence();
    }
}
```

Once registered, this function is now available across the whole platform, and can be used in multiple contexts.
For instance, it can be found under `Help | Functions`, and executed right from there: 

![](complement.png)

Since this function is tagged as `converter` and declares semantic type of the input parameter, it 
becomes a converter. It will be offered to convert columns with the semantic type "nucleotide". Note that
the semantic type detector function can be registered in the same package.  

![](complement-converter.png)

While the function is available in the platform's UI, the corresponding JavaScript file and
its dependencies do not get loaded in the browser until a function is used for the first time. This allows the
platform to operate on thousands of functions covering different domains without compromising the performance.    

Each package has a version number. All objects inside the package are being deployed according to package version. 
That means that if version changes, there will be an independent instance of each package asset.  
Multiple versions of each package can be deployed at one moment, and administrator can switch between them. All users 
will see only objects that belong to current package version.
There is a special "debug" version can be deployed for each package. If developer deploys a debug version, it becomes 
active for current package until developer deletes it or changes his developer key. In this case developer can see objects 
from his version of package, and changes don't affect other users package representation.  


## Buildable Packages

To get started, you need NodeJS and NPM installed. Also, install WebPack to be able to build package locally
and debug it using WebPack Dev Server.

To get a package template:
1. Create an empty folder
2. Run `npm install datagrok-tools -g`. It's a tool to help you to build and deploy your package.
3. Run `datagrok-init`. Enter package name, remote Datagrok server URI (including /api) and your developer key.
This will drop to yor package all necessary files. You'll be able to change all of this in the future.
4. Run `npm install`. This will install packages according to package.json dependencies list. For the first time, 
there are only the `datagrok-api` package, and you are free to add new dependencies as usual.

Take a look at the `src/package.js` file:
```js
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name test
//input: string s
export function test(s) {
    grok.shell.info(_package.webRoot);
}
```
Note, that there are already Datagrok API modules imported. They are also set as external modules, so Webpack 
doesn't include them to the output.
You are completely free to include any other libraries or packages, as they are all will be build in the single bundle file.
To debug package run `npm run upload-debug`. This will upload package to the server in debug mode, just for your account.
To deploy package for all users run `npm run upload-deploy`. Also, package could be deployed from git, or another source. 

To change Datagrok server - edit scripts section in the `package.json` file and add new dev key to the `upload.keys.json`. 

## Development

A JavaScript package runs inside the Datagrok platform, which makes the development and
debugging experience different compared to the more traditional web applications. Essentially, the packages
are developed locally, but the platform runs remotely. To enable the best possible experience for
developers, we established a workflow where the package is uploaded to the remote server at startup, and
then gets served from the server. By associating local JavaScript files with the remote sources in your favorite IDE, 
it is possible to hide the complexity of that scenario. For instance, you can set breakpoints, do 
step-by-step execution and generally debug the program in the regular way. Of course, you can always
use the debugger that comes with the browser.

To develop Grok packages, we recommend to start with downloading one of the existing package templates.
Then, set up your IDE in such a way that when starting a project, it would [publish](#publishing) the package,
and then start the platform. Here is how to do it in WebStorm:

TODO: insert gif   

Packages deployed in the development mode are visible only to the authors. This ensures that multiple
people can simultaneously work on the same package.    

### General notes on package development

Our approach to extending the system is providing one canonical, convenient way for 
developers to achieve the task, at the same time exposing enough extension points
to enable deep customization. Same concepts apply to JavaScript development. 
We do not impose any requirements on the UI frameworks or technologies used for 
the JavaScript plugins, although we encourage developers to keep it simple.

To simplify development, Datagrok provides an "Inspector" tool that lets developers peek
under the hood of the platform. Use it for understanding which events get fired and when,
how views and viewers are serialized, what is getting stored locally, what widgets are
currently registered with the system, etc. 

### Environments

In order to isolate packages being debugged from the production instance, we recommend running 
them against the "dev" instance, if possible.

## Publishing

When a package gets published, a "published package" entity gets created. It is associated with the 
package, and has additional metadata (such as published date). Typically, only one version of a package
is visible to a user. Administrators can manage published packages and decide which versions
should be used. It is possible to roll back to an old version, or assign a particular version to
a particular group of users.  

To publish a package, use the "publish" tool included in the SDK. Alternatively, drop a zip 
file with the package content into the platform. The tool supports two deployment modes:

* Development - the package will only be accessible by the developer
  ```upload_package -p <package_path> -k <dev_key> -r <api_endpoint>```
* Production - available to everyone who has the necessary privilege
  ```upload_package -p <package_path> -k <dev_key> -r <api_endpoint> deploy```

[Video guide to package publishing](https://www.youtube.com/watch?v=PDcXLMsu6UM)

Package template contains two cmd files, "upload-debug.cmd" and "upload-deploy" that have all
parameters configured already (including your dev key). Dev key can also be retrieved by opening your user profile
and clicking on `Developer key`. Administrators can manage existing keys and grant or revoke developer
privileges.

To see available packages, click on `Admin | Packages`, or follow [this link](https://public.datagrok.ai/packages)
from outside the platform. Just as pretty much anything in the platform, packages are subject to
[privileges](../govern/security.md#privileges), which allows a fine-grained control over the
platform's content.  

## Applications

Applications are [functions](../overview/functions/function.md) tagged with the `#app` tag. A package might contain zero, one, or more apps. 
To open application launcher, click on `Admin | Apps`, or follow [this link](https://public.datagrok.ai/apps)
from outside the platform.  

To launch a particular app automatically, open the following URL: `https://public.datagrok.ai/apps/<APP_NAME>`

## Documentation

According to [this study](http://sigdoc.acm.org/wp-content/uploads/2019/01/CDQ18002_Meng_Steinhardt_Schubert.pdf),
in terms of the strategies used for understanding API documentation, different developers fall into few groups 
(systematic, opportunistic and pragmatic). These findings are consistent with our experience. 
For Datagrok's documentation, we have established an approach that enables developers from either of the
above-mentioned groups to be productive. 

* [Sample browser](https://public.datagrok.ai/js) (`Tools | Scripting | JavaScript`) is an interactive tool
  for browsing, editing, and running JavaScript samples that come with the platform. Samples are grouped by
  domain, such as data manipulation, visualization, or cheminformatics. They are short, clean examples
  of the working code using [Grok API](js-api.md) that can be copy-and-pasted into the existing solution.
  The samples are also cross-linked with the [help](https://datagrok.ai/help) system. 
* [Grok API](js-api.md) provides complete control over the platform. 
  [JSDoc documentation](https://datagrok.ai/help/develop/api/index.html) is available.
* [Platform help](https://datagrok.ai/help/) explains the functionality from the user's point of view. Where
  appropriate, it is hyper-linked to samples and demo projects. In the near future, we plan to turn it
  into the community wiki, where users would be contributing to the content. The same web pages
  are used as an interactive help within the platform (you see help on the currently selected object).
  
Additionally, there are few ways to connect with fellow developers:    
* [Datagrok forums](https://public.datagrok.ai/forums) for in-depth, topic-centered discussions. 
  It is possible to include "live" content in messages, such as viewers, project, users, etc. 
  The UX is still a "work in progress", but it's maturing quickly.   
* [Slack space](https://datagrok.slack.com) 

See also: 
* [Grok API](js-api.md)
* [How Developers Use API Documentation: An Observation Study](http://sigdoc.acm.org/wp-content/uploads/2019/01/CDQ18002_Meng_Steinhardt_Schubert.pdf)

## Roadmap

* Developing shareable components and widgets

### Videos

<iframe width="560" height="315" src="https://www.youtube.com/embed/PDcXLMsu6UM" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

See also

* [Github repository](https://github.com/datagrok-ai/public/tree/master/packages)