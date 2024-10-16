# Datagrok-tools changelog


## 4.13.31 (2024-10-15)

### Features

* Grok test without repository fixes


## 4.13.30 (2024-10-11)

### Features

* Grok Check fixes

## 4.13.29 (2024-10-09)

### Features

* Help fixes


## 4.13.28 (2024-10-09)

### Features

* Grok Check update(added checks on caching)

## 4.13.27 (2024-10-02)

### Features

* Publish all added

## 4.13.26 (2024-10-02)

### Features

* Create command fixes

## 4.13.25 (2024-09-30)

### Features

* Test command package publishing fix 

## 4.13.24 (2024-09-26)

### Features

* Test command fixes update

## 4.13.23 (2024-09-23)

### Features

* Test all bug fixes

## 4.13.22 (2024-09-20)

### Features

* Implemented mechanism which allows to test few packages by one function

## 4.13.21 (2024-09-17)

### Features

* Grok check updated for release candidate versions

## 4.13.20 (2024-08-28)

### Features

* Grok link command fixes 

## 4.13.19 (2024-08-27)

### Features

* Package ts data insertation fixes

## 4.13.18 (2024-08-27)

### Features

* Creation fixes

## 4.13.17 (2024-08-27)

### Features

* Added decorators for viewer functions

## 4.13.16 (2024-08-23)

### Features

* Grok link fixes

## 4.13.15 (2024-08-21)

### Features

* Grok link updated to link all dependent libs 

## 4.13.14 (2024-08-09)

### Features

* Timeout update for test runs up to 1 hour

## 4.13.13 (2024-08-09)

### Features

* Added stress test flag

## 4.13.12 (2024-08-05)

### Features

* Puppeteer-screen-recorder version fixed to 3.0.3

## 4.13.11 (2024-07-30)

### Features

* Added commit label to context panel for 'manually' published versions

## 4.13.10 (2024-07-29)

### Features

* Added ability to select test by args variable 

## 4.13.9 (2024-07-29)

### Features

* Added ability to select category by args variable 

## 4.13.8 (2024-07-18)

### Bug Fixes

* Fixed path chacks for npmignore

## 4.13.6 (2024-07-15)

### Features

* Splitted warnings and errors in grok check command

## 4.13.4 (2024-06-25)

### Bug Fixes

* Check source maps fixes 

## 4.13.3 (2024-06-21)

### Features

* Added check on source maps in grock check command

## 4.13.1 (2024-06-04)

### Features

* Added argument `package`, which allows you to choose a package for the `test` command


## 4.12.13 (2023-08-17)

### Bug Fixes

* Datagrok API version check fix
* Check Changelog fix
* Sync --csv and --verbose flags

## 4.12.12 (2023-08-07)

### Features

* Check for datagrok-api dependency

### Bug Fixes

* Latest package version in CHANGELOG check fix

## 4.12.11 (2023-08-04)

### Features

* GROK-13643 Check improvements:
  * There is no beta property in package.json
  * No datagrok-tools in dependencies (or latest version)
  * Latest version from package.json is in CHANGELOG (warning)
  * Ignore CHANGELOG checks for service packages ("servicePackage" property in package.json)
  * Change supported h2 formats (1.7.9 (2023-07-24) and 1.7.9 (WIP))
  * For packages < 1.0.0 exit with exit code 0, and only show warnings. And for packages >= 1.0.0, exit with a non-zero code (only for check command)
  * If an invalid flag/command is specified, output the help and exit with exit code 1

## 4.12.10 (2023-08-01)

### Features

* Video recording enhancements

### Bug Fixes

* FuncSignatures check fix

## 4.12.7 (2023-07-28)

### Features

* Tools: Changelog h2 new format

## 4.12.4 (2023-07-24)

### Features

* GROK-13573 Tools: simplify output (add --verbose flag)
* GROK-13573 Tools: checks for changelog
