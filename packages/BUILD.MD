# Build instructions

This repository contains 3 types of packages:
1. JS API itself: `/js-api`
2. Libraries: `/libraries`
3. Regular Datagrok packages: `/packages`

# Building JS-API and libraries

To build JS-SPI simply run `npm install` and `npm rub build` inside package folder. It builds Typescript code to Javascript.

# Working with package dependencies

Each package or library must have a certain list of dependencies. It can depend on js-api or library.
Usually, you must specify minimal library version package can work with, otherwise, you can't build a package because of lack of code it needs.
If package uses features from js-api that implemented in the new version, you must specify correct version constraint, i.e. `>=0.108.0`
If package uses features from libraries, you must specify correct version constraint for library.
Run `npm install` inside a package folder, and Node downloads all dependencies to `node_modules` folder of the package, including Datagrok and third-party packages.

Feel free to change package code and publish it anywhere with `grok publish` for debugging and testing.

# Publishing packages and libraries

Publishing the package or the library to NPM is the right way to make you changes available for customers and other developers.

To publish something you need to bump package or library version and push this change to the repository.
Please, ensure that GitHub action completed successfully. 

If you add something to library, it is your responsibility to publish it.
If you add something to js-api, please, note that others can use it only after publishing of a new platform version.

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
