# Build instructions

This repository contains 3 types of packages:
1. Datagrok JS API: `/js-api`
2. Shared libraries: `/libraries`
3. Datagrok packages: `/packages`

# Building JS-API and libraries

To build JS-API, run `npm install` and `npm pub build` inside the package folder. It builds Typescript code to Javascript.

# Working with package dependencies

Each package or library must have a certain list of dependencies. It can depend on js-api or on a library.
Make sure to  specify minimal library version a package can work with.
If a package uses features from js-api that are implemented in the new version, you must specify the version constraint, i.e. `>=0.108.0`
If a package uses features from libraries, you must specify the correct version constraint for library.

To install package dependencies, run `npm install` from the package folder. 

To publish a package, run `grok publish` from the package folder. See [this guide](../help/develop/develop.md) for 
details.

# Publishing packages and libraries

To make you changes available to others, publish a package to NPM. 

To publish a package or a library, you need to bump package version, and push this change to the repository.
This will trigger the automatic GitHub action. Make sure that GitHub action completed successfully before
trying to use the new package. 

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

There is a way to skip NPM publish step for local debugging. It might be useful if you modify the code in many libraries simultaneously and want to debug it. 
If you need to use local copy of a library or js-api as a dependency of your package instead of NPM version, you need:
```
```
