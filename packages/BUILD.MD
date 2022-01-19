# Build instructions

This repository contains 3 types of packages:
1. JS API itself: `/js-api`
2. Libraries: `/libraries`
3. Regular Datagrok packages: `/packages`

# Building JS-API and libraries

To build JS-SPI simply run `npm rub build` inside package folder. It builds Typescript code to Javascript.

# Package dependencies

Each package must have correct dependencies in `package.json`. 
If package uses features from js-api that implemented in the new version, you must specify correct version constraint, i.e. `>=0.108.0`
If package uses features from libraries, you must specify correct version constraint for library.

Feel free to change package code and publish it anywhere with `grok publish`. 

# NPM scopes

Datagrok has 2 scopes

## @datagrok

Scope for public Datagrok packages. All Datagrok instances around the world publish packages from this scope. 
All packages from this GitHub repository automatically go to NPM repository right after you change package version.
If package builds successfully and passes tests, it becomes available in this NPM scope.

## @datagrok-libraries

Repository for libraries. 
Just like the package, if you change the version of library and if it passes tests, library becomes available in NPM scope

