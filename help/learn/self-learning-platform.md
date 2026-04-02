---
title: "Self-learning platform"
---

One of Datagrok's unique features is that most data operations happen within the
platform. Beyond the convenience of a single tool with centralized access management,
this lets the system learn from observed user behavior.

Think about Netflix's movie recommendation engine, but instead of dealing with just two entities (
users and movies) and one relation (user's score for the movie) we have a much more complex case. We got dozen of
different [entity types](../datagrok/concepts/objects.md)
(such as [query](../access/access.md#data-query), [viewer](../visualize/viewers/viewers.md), etc), connected with different relations
(such as '[query](../access/access.md#data-query) `ran_by` [user](../govern/access-control/users-and-groups#users)') and restricted by different constraints.

When enabled, the self-learning component uses AI to spot usage patterns and suggest
relevant actions. For example, it might suggest a visualization for your dataset,
apply a [prediction model](learn.md) trained by a colleague, or create a derived
column like BMI from weight and height columns. Suggestions draw from all users'
activity, not just yours, helping spread organizational knowledge across departments
and time zones.

See also:

* [Predictive modeling](learn.md)
* [Data queries](../access/access.md#data-query)
* [Spaces](../datagrok/concepts/project/space.md)
