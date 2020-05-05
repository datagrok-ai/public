<!-- TITLE: User -->
<!-- SUBTITLE: -->

# User

Represents the identity of a person, and is used for [group](group.md) and role management. 

## Filtering

Within the users view, you can use these fields to filter users with [smart search](../overview/smart-search.md):

| Field       | Description                                        |
|-------------|----------------------------------------------------|
| id          |                                                    |
| name        |                                                    |
| login       |                                                    |
| email       |                                                    |
| createdOn   |                                                    |
| updatedOn   |                                                    | 
| group       | [Group](group.md) object: User personal security group |
| author      | [User](user.md) object                                 |
| starredBy   | [User](user.md) object                                 |
| commentedBy | [User](user.md) object                                 |


## Profile

User profile view contains summary information about a particular user. To open it, click on your user's icon 
in the right top corner in order to open your profile page. To open a profile page of another user, 
right-click on the name, and select "Details". User profile contains several panes.

**Summary** pane contains important messages retrieved from the [audit](../govern/audit.md).
Additionally, there are links for  changing your password, 
obtaining [developer key](../develop/develop.md#publishing), API key, and logging out of the platform.  

**Chats** pane contains chats that you have participated in.

**Activity** pane lists user actions recorded in the platform, with the most recent actions on top. 
It uses [audit](../govern/audit.md) as a data source. Note that the activity is interactive, i.e. 
by clicking on highlighted entities their properties will appear in the [property panel](../features/property-panel.md). 

**Favorites** contains your favorite objects. To add something to favorites, click on it, and 
then click on the "star" icon on the right in the property panel. Not every object can be added to favorites,
currently the list is limited to connections, queries, projects, and functions.  

**Projects**, **Connections**, **Queries**, **Models**, **Scripts**, **Projects**, and **Notebooks**
contains corresponding objects that were created or used by you.


See also:

* [Group](group.md)
* [Audit](audit.md)
* [Authentication](authentication.md)
* [Security](security.md)
