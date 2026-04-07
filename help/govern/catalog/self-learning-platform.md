---
title: "Self-learning platform"
---

One of Datagrok's unique features is the ability to perform most data operations from within the
platform. Besides the convenience of keeping everything in one place, having centralized access management,
and reducing the number of tools, this unlocks the capability for the system to learn based on observed user
behavior.

Think about Netflix's movie recommendation engine, but instead of dealing with just two entities (
users and movies) and one relation (user's score for the movie), you have a much more complex case. There are dozens of
different [entity types](../../datagrok/concepts/objects.md)
(such as [query](../../access/access.md#data-query), [viewer](../../visualize/viewers/viewers.md), etc.), connected with different relations
(such as '[query](../../access/access.md#data-query) `ran_by` [user](../access-control/users-and-groups.md#users)') and restricted by different constraints.

When enabled, the self-learning component uses various AI techniques to spot patterns in usage and provide users with
actionable insights. These might include suggestions to visualize the currently open dataset in a specific way, predict
properties based on a [prediction model](../../learn/learn.md) trained by someone else, create a derived column (such as
BMI if your dataset contains weight and height), or many other actions. The platform gives you
suggestions based not solely on your activity but on the activity of other users as well. This facilitates
spreading your organization's knowledge across different departments and time zones.

See also:

* [Predictive modeling](../../learn/learn.md)
* [Data queries](../../access/access.md#data-query)
* [Dashboards](../../datagrok/concepts/project/dashboard.md)
* [Spaces](../../datagrok/concepts/project/space.md)
