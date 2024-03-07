# Grok Connect changelog

## 2.1.14 (2024-02-21) (WIP)

* Increase WebSocket max size of incoming string message

## 2.1.13 (2024-01-19)

### Bug fixes

* Fix empty DataFrame with no columns when no results were returned from the query
* Fix MSSQL credentials expose in logs
* Fix TableQuery. Creates a query with aliases for column names
* Fix blob, clob and nclob display for Oracle
* Fix deadlock when several drivers have been loaded in parallel

## 2.1.12

* Make logging async
* Improve queries logging and debugging

## 2.1.10 - 2.1.11

* Fix categorize()
* Fix NullPointer when update query executed
* Fix Snowflake connection building
* Fix bit type support in columns

## 2.1.9

* Update of Neptune driver version

## 2.1.8 

This release focuses on fixing bugs

### Bug fixes

* Fix connection form of Athena Provider
* Fix exposing sensitive data in logs

## 2.1.7 (2023-07-28)

This release focuses on fixing bugs

### Features

* Add possibility to set log level in command line args when grok connect is starting

### Bug fixes

* Fix connection form of Impala Provider
* Fix datetime parameter set in UTC
* Remove connection pool for Neptune due to memory leak

## 2.1.6 (2023-07-24)

This release focuses on fixing bugs

### Bug fixes

* Fix bug with socket data sending with the size of 100 when fetch size change is not supported
