# Build instructions

This repository contains 3 types of packages:
1. Datagrok JS API: `/js-api`
2. Shared libraries: `/libraries`
3. Datagrok packages: `/packages`

# Building libraries and packages

To build a library or a package run `npm install` and `npm run build` inside the package folder. 
It builds Typescript code to Javascript bundle.

# Building js-api (datagrok-api package)

There is no need to build it separately except for local debugging (see [linking](#linking)). CI system builds it with the platform. 
If you still want to test your code with custom js-api package, you can [inject it](#api-injection)

# Working with package dependencies

Each package or library must have a certain list of dependencies. It can depend on js-api or on a library.
Make sure to specify minimal library (or js-api) version a package can work with.
If a package uses features from js-api that are implemented in the new version, you must specify the version constraint, i.e. `^0.108.0`
If a package uses features from libraries, you must specify the correct version constraint for library.

To install package dependencies, run `npm install` from the package folder. 

# Deploying packages to Datagrok server

To publish a package to Datagrok server, run `grok publish` from the package folder. See [this guide](../help/develop/develop.md) for
details.
Avoid using `--rebuild` key, always try to build package locally using `npm run build`
 
# Publishing packages and libraries to NPM

To make you changes available to others, publish a package to NPM. 

To publish a package or a library, you need to bump package version, and push this change to the repository.
This will trigger the automatic GitHub action. Make sure that GitHub action completed successfully before
trying to use the new package. Package version must comply with Semantic Versioning.

It is the package author's responsibility to publish the package/library when the time is right.

If you add something to js-api, note that the changes will be available with the next core release.

# NPM scopes

Datagrok has 2 scopes

## @datagrok

Scope for public Datagrok packages. All Datagrok instances around the world publish packages from this scope. 
All packages from this GitHub repository automatically go to NPM repository right after you change package version.
If package builds successfully and passes tests, it becomes available in this NPM scope.

## @datagrok-libraries

Repository for libraries. 
Just like the package, if you change the version of library and if it passes tests, library becomes available in NPM scope

# Linking

[NPM link](https://docs.npmjs.com/cli/v8/commands/npm-link) is a way to skip NPM publish step for local debugging. It might be useful if you modify the code in many libraries simultaneously and want to debug it. 
If you need to use local copy of a library or js-api as a dependency of your package instead of NPM version, you need:
1. Run `npm link` in library (or js-api) folder
2. Run `npm link LIBRARY_NAME` in package folder. For example: `npm link @datagrok-libraries/utils`
That makes NPM use local library instead of NPM version. To reset link, simply run `npm install` in package folder.
Usually, there is a script in a package or in a library, called `npm run link-all` to link all local dependencies. 
To link all the libraries run `libraries-link-local.cmd` script. Refer to script contents to understand what's going on.

> Suppose, you need to debug package `Chem` with local versions of `js-api` and `utils` library
> 1. Run `libraries-link-local.cmd`.
> 2. Run `npm run link-all` in `Chem` package folder.
> 3. Run `npm run build-chem-local` to build all libraries and package itself.
> To revert changes and make package use NPM libraries simply run `npm install`.
> Don't forget to bump library version, package version and package dependencies properly before commit.

You can publish built package to any Datagrok server for testing and debugging, but if js-api version, the package was built with, differs from platform version - package might be unstable.
> <span style="color:red">!</span> If you want to link multiple dependencies - link them in a single command.

> <span style="color:red">!</span> When you link package to local js-api, don't forget to link all libraries to local js-api too

# API injection

## Platform injection 

There is a way to inject your js-api inside local Docker container. Use [Docker Volumes](https://docs.docker.com/storage/volumes/) to map folder with pre-built js-api to your local folder
1. Run `npm run build` inside `js-api` folder. By default, it saves a bundle to `../../xamgle/web/js/api` folder. Remember this path.
2. Map this folder to Datagrok Docker:
 - Open the `localhost.docker-compose.yaml` file (in `public\docker` folder).
 - Replace 
```
volumes:
      - datagrok_data:/home/grok/data
      - datagrok_cfg:/home/grok/cfg
``` 
section by the following code
```
    volumes:
      - datagrok_data:/home/grok/data
      - datagrok_cfg:/home/grok/cfg
      - type: bind
        source: %FULL_PATH_TO_BUILD_JS_API%
        target: /home/grok/datagrok/html.build/web/js/api
```
3. Restart docker containers as it is described in local Docker [instructions](https://datagrok.ai/help/develop/admin/docker-compose)
4. Clear browser cache, refresh browser page. 

## Package injection

Remove `datagrok-api` from externals section of your `webpack.config.js` file:
```
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
```
This will include js-api inside package bundle instead of using platform code. Don't forget to revert this change before publishing.
