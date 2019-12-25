<!-- TITLE: Audit -->
<!-- SUBTITLE: -->

# Audit

Audit system is intended to store all users activity for future analysis. 

All of the changes made to the different [entities](../entities/entities.md), is tracked and can be analyzed using activity sections.
All events are joined to user session. Changes, made to the entities, are connected with corresponding entities.
Audit records can be posted both from client and server process.

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