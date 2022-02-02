<!-- TITLE: Self-learning platform -->
<!-- SUBTITLE: -->

# Self-learning platform

One of the Datagrok's unique features is the ability to perform most of the operations on the data from within the
platform. Besides the obvious convenience of keeping everything in one place, having a centralized access management,
and reducing the number of tools, this unlocks the capability for the system to learn, based on the observed user's
behavior.

Think about Netflix's movie recommendation engine, but instead of dealing with just two entities (
users and movies) and one relation (user's score for the movie) we have a much more complex case. We got dozen of
different [entity types](../overview/objects.md)
(such as [query](../access/data-query.md), [viewer](../visualize/viewers.md), etc), connected with different relations
(such as '[query](../access/data-query.md) `ran_by` [user](../govern/user.md)') and restricted by different constraints.

When enabled, the self-learning component uses various AI techniques to spot patterns in usage, and provide users with
actionable insights. These might include suggestions to visualize currently open dataset in a specific way, predict
properties based on [prediction model](predictive-modeling.md) trained by someone else, create a derived column (such as
BMI in case your dataset contains weight and height), or many other actions. Note that the platform would give you
suggestions that are based not solely on your activity, but on the activity of other users as well. This facilitates
spreading organization's knowledge across different departments and time zones.

See also:

* [Predictive modeling](predictive-modeling.md)
* [Data queries](../access/data-query.md)
* [Projects](../overview/project.md)
