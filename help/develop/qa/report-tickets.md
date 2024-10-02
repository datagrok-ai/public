---
title: "Report bugs and feature requests"
---

This guide provides a set of rules for reporting issues and managing
tickets in Datagrok. We use two issue tracking systems:

* [JIRA](https://reddata.atlassian.net/)
for internal issues integrated with BitBucket
* [GitHub Tracker](https://github.com/datagrok-ai/public/issues)
for all externally reported issues.

We recommend reporting issues in JIRA if you have access to it, but if you
don't, use GitHub Tracker. The rules of the bug report format adopted by the
company do not differ in both issue tracking systems.

## Reporting an issue

A good ticket should have a type, title (summary), description, assignee,
priority, and due date. Some tickets should be included in the projects or
epics. Here we provide recommendations for each item:

* **Issue type**. In GitHub, we use two labels to type an issue: _bug_ and
  _enhancement_. In Jira, we use similar issue types: _bug_, _improvement_
  (enhancement of existing functionality or UI), and _new feature_ (requesting
  new capability or software feature).
* **Title (Summary)**. Start the title with the platform module or package name,
  then name unexpected behavior for a bug and a new feature name or improvement SWsummary for an enhancement. 

|Issue type|Formula|Example|
|----------|-------|-------|
|Bug report| [Platform module] or [Package name] : [Unexpected behavior]| Line chart: Y-axis label is broken|
|New feature or Improvement| [Platform module] or [Package name] : [New feature name] or [Improvement summary]| Bar Chart: context menu harmonization|

* **Description** is important for testing and writing release notes. For an
enhancement, write the description of the new functionality. A bug description
should include the following items:
  * **Instance and its version**. If the bug concerns a package functionality,
  mention the package version.
  * **Data on which the bug is reproduced (optional)**. If a bug occurs on a
  particular dataset only, include it in the bug description.
  * **Steps to reproduce**. A step-by-step description of your actions causing
  the bug (where to click, what to enter into inputs, which checkboxes to
  mark, etc.)
  * **Expected Result**. Describe how the platform is supposed to work—the ideal
  result you want but don't get.
  * **Actual result**. The result that you receive when using the platform,
  which you don't expect to receive.
  * **Code examples (optional)**. If the bug concerns working with the platform
  API, scripts, queries, or something where you use a code, give an example.
  * **Priority**. Set the priority for the bug in Jira or label the issue in
    GitHub. Use your common sense above all rules described below.
    * _Highest_. It is a bug in the main functionality of the platform. Such as
    Chem Descriptors, Detectors, and so on.
    * _High_. An important request from the customer, something that blocks
    others' work, or a task with the exact due date. Set the due date for such
    tasks even if it was not initially specified by the issue's reporter. In
    most cases, such requests should be resolved in a week.
    * _Medium_. The default priority for the task.
    * _Low_. A good enhancement for the platform which does not have a big
      impact on the users.
  * **Label (optional)**. Use labels to display the ticket on the relevant board
    or project.
  * **Attachment (optional)**. You are welcome to attach screenshots, GIFs,
  recorded videos, etc. to complete the issue description.
* **Due date**. Set the due date according to the priority, client's request, or
corresponding workstream date. Generally, medium-priority issues should be
done within two weeks, high-priority within a week, and the highest-priority
tickets within a few days.
  >Note: In GitHub, you are allowed to set the due date only for a ticket within
  >a project. Thus for stand alone tickets in GitHub, adjust the due date with
  >labels of priority (high, low).
* **Assignee**. Assign Alexander Paramonov or Denis Kryvchuk if the issue
relates to the platform core. To understand the person responsible for a
particular package, look inside the package.json file, which must exist for
each package. Inside this file, there is an "author" field from which you can
find out whom to assign.
* **Epic link (optional)**. Use it to include a ticket in the epic.

## Managing tickets

To track the actual state of the platform, it's important to keep tickets' status up to date. Also, use status to interact with other participants in the process:

* Set **In Progress** status when you start working on the issue.
* Use **Testing** or **In Review** status when you or your college is testing or reviewing your work on the issue.
* Set the status to **Ready for QA** if you think your work is done and it's time for the QA engineer to test it.
If there are no errors or comments, the QA engineer sets the Done status for the ticket.

>Note: Write meaningful commits, and don't forget to mention the issue number. For details, see the [commit message policy](../dev-process/git-policy.mdx).
