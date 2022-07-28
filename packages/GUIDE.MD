# Package ownership

Each package must have a responsible person.

This person DOES:
- maintain the package
- keep package tests green 
- gather issues and tickets about the package
- refactor the package
- check and test package with all new Datagrok releases

This person DOES NOT:
- write all the package code
- fix all the package problems

This person IS:
- always aware of package status
- an entry point for customers questions and issues
- in charge of test coverage
- in charge of GitHub Action status and publication status

This person IS NOT:
- the only package developer

## Package files and metadata
### Naming and Versioning

Each package must be properly named.

* Use kebab-case for package name (`power-pack`) 
* Use PascalCase for package folder name (`PowerPack`).
* Use kebab-case for file names (`compare-columns.ts`).
* Name must contain scope `@datagrok` either `@datagrok-libraries`.

Examples: 
```
@datagrok/power-pack
@datagrok/clinical-case
@datagrok/chem
@datagrok-libraries/utils

MySuper-package !!!INCORRECT
```

Each package must contain a version tag in `package.json` file. 
Version must follow [Semantic Versioning](https://semver.org/) notation.
Use `0.0.1` version as a starting point.
When package is ready to be published widely, set 1.0.0 version.

### Documentation

Each package must contain a short description in `package.json` and a `README.MD` file.
Each package must contain `author` section in `package.json`:
```
"author": {
  "name" : "Alexander Paramonov",
  "email" : "aparamonov@datagrok.ai"
}
```

See also: [NPM package.json](https://docs.npmjs.com/cli/v8/configuring-npm/package-json)

### Repository

`package.json` must contain:
```
  "repository": {
    "type": "git",
    "url": "https://github.com/datagrok-ai/public.git",
    "directory": "packages/CustomML"
  }
  ```

### Dependencies

Dependencies section with `datagrok-api` dependency is mandatory. 
Refer to [Semantic Versioning](https://semver.org/) concept to properly set `datagrok-api` version constraint.
All "compile-time" packages must be in `devDependencies` section.:
```
	"devDependencies": {
		"@types/jest": "^27.0.0",
		"jest": "^27.0.0",
		"jest-html-reporter": "^3.6.0",
		"path": "^0.12.7",
		"puppeteer": "^13.7.0",
		"ts-jest": "^27.0.0",
		"ts-loader": "^9.2.6",
		"typescript": "^4.4.4",
		"webpack": "^5.59.1",
		"webpack-cli": "^4.9.1"
	},
```

### Package Sources
  
Avoid using `sources` section in `package.json`. Always prefer to link external package as a dependency and import as a module.
If it's impossible, first try to use `common/` package. Place external packages files outside `/src` folder, `/vendor` is a good place.

### Directory structure

Basic package structure should look like that
```
  css/
    ...
  files/ 
    ...
  scripts/
    ... 
  src/
    __jest__/
      remote.test.ts
    tests/
      ...
    package.ts
    package-test.ts
  vendor/
    ...
  package.json
  package.png
  .npmignore
  jest.config.js
  README.MD
  tsconfig.json
  webpack.config.js
```

## Tests

Each package must have auto-tests. We use custom Test Framework that is similar to Jest, and it supports Datagrok JS-API.
One single Jest test must be set up to enable package testing.
Please, refer to `Chem` package, it has Jest properly configured and all tests written properly.

You must always keep tests green, you can't publish a package without all tests passing.

## Package content
 
Package doesn't just exist. 

It is used to solve the customer problem and not to add another one. Always keep that in mind. 

Be sure that the package:
- has an autostart function only if it must have an autostart function.
- is small
- is fast
- doesn't break the platform
- doesn't break other packages
- doesn't conflict with other packages
- looks cool and fancy

Again, the responsible person is not always a person who fixes all the problems, but the person, who is always aware and controls the process.

## Code style


## Troubleshooting