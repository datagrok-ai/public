<!-- TITLE: Audit -->
<!-- SUBTITLE: -->

# Audit

Audit system is intended to store all users activity for future analysis.

All of the changes made to the different [entities](../datagrok/objects.md), is tracked and can be analyzed using
activity sections. All events are joined to user session. Changes, made to the entities, are connected with
corresponding entities. Audit records can be posted both from client and server process.

Each audit event can be disabled remotely from server for certain users group.

All entities have activity section in property panel with all corresponding audit records.

All audit records have fixed types:

* query-created
* query-edited
* query-deleted
* query-start
* query-published
* query-transformations-edited
* connection-created
* connection-edited
* connection-deleted
* connection-published
* job-created
* job-edited
* job-deleted
* job-transformations-edited
* job-start
* script-created
* script-edited
* script-deleted
* script-start
* script-published
* predictive-model-created
* predictive-model-edited
* predictive-model-deleted
* predictive-model-start
* predictive-model-published
* action-start
* notebook-created
* notebook-edited
* notebook-deleted
* notebook-start
* notebook-published
* project-created
* project-edited
* project-deleted
* project-opened
* entity-shared
* entity-shared-silent
* table-produced
* user-invited
* comment-posted
* dialog-ok
* main-menu-item-click
* error
* tutorial-completed
* viewer-rendered
* log
* package-tested

## JavaScript API

You can create a Logger object to put Audit Records to Datagrok server:

```javascript
// Create logger
let logger = new DG.Logger();
// You can pass a callback that sets predefined parameters
let logger = new DG.Logger((m) => m.params['persistent']= 'value');
```

Use `logger.log()` method to save records:

```javascript
// Default type is "log"
logger.log('HELLO WORLD', {test: 'value', 'foo': 'bar'});
// But you can specify another
logger.log('HELLO WORLD', {test: 'value', 'foo': 'bar'}, 'my_log');
```

All errors are automatically logged to Datagrok, but you can log store them explicitly:

```javascript
// Specify Stack Trace to be able to group errors
logger.log('Error!', {'stackTrace': '...'}, 'error');
```

## Amazon cloud watch export

Configure automated export Amazon Cloud Watch Log in `Settings -> Log`.

Specify settings:

### Cloud watch log group

Fill with AWS CWL log group name

### Events to cloud watch streams map

Choose which audit types go to log streams. For example `log => datagrok_log`
, `error => datagrok_errors`, etc.

### Cloud watch access key, cloud watch secret key, cloud watch region

Amazon Access and secret keys and region
