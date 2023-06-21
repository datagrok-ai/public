---
title: "Content types"
---

We follow four different content types:

* [Tutorials](#tutorials)
* [How-To guides](#how-to-guides)
* [Technical reference material](#reference-material)
* [Topics](#topics)

These content types represent different documentation functions and require a distinct writing mode. Users need all four
kinds of documentation at different times, in different circumstances.

|         | Tutorials                  | How-To guides            | Technical reference                                               | Topics                                                                                                            |
|:--------|:---------------------------|:-------------------------|:------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------|
| Purpose | Learn                      | Achieve a goal           | Get information                                                   | Understand                                                                                                        |
| Goal    | Allow users to get started | Solve a specific problem | Describe technical aspects of how things work and how to use them | Describe key concepts and topics at a fairly high level and provide useful background information and explanation |
| Form    | A walk-through             | A series of steps        | Dry description                                                   | Explanation                                                                                                       |

Some of the Datagrok documentation content doesn't classify as one of the four primary content types:

* _Release notes_
* _FAQs_ are there to answer questions that go beyond the scope of tutorials and How-To guides.
* _Internals_ is where you find information about how to contribute code, report bugs, request features, or write
  documentation.

## Tutorials

Tutorials walk the user through a series of steps to create something. Tutorials explain the nature of the problem so
that the user understands what they are trying to achieve. Like How-To guides, tutorials are action-oriented. Don't
explain why a user needs to do something. Tell them how to do it instead.

In general, tutorials are most useful when users are just getting started. For example:

* Tutorials about Datagrok UI.
* Tutorials about a complex workflow with multiple steps and sub-steps.
* Tutorials about a Datagrok feature, capability, or package (such as _Cheiminformatics tutorial_ or _Exploratory data
  analysis tutorial_).

Try to limit tutorials to 10 steps or 30 minutes in length. More than that may indicate you're trying to describe a use
case that is too complex for beginners or combining too many actions in one tutorial. The more complex and
time-consuming a tutorial is, the less likely the user will finish it.

Tutorials should be in this format:

1. Start by telling the user what the tutorial does, the estimated completion time, and the expected _outcome_:
    * _By the end of this tutorial, you will be able to_ followed by a phrase that starts with a verb and clearly
      describes the outcome of what the user will accomplish.
    * _Estimated completion time: 30 minutes._
2. Provide prerequisites:
    * _Before you begin, make sure you meet these prerequisites:_
    * _A conceptual understanding about x._
    * _Access to the y application._
3. List everything the user will do.
4. Next, tell users how to complete each task. Group related tasks into sections. Set clear expectations for each step
   or task.

   Use this template:

    * _Step 1._ Write the command here. Start with the verb.
    * Add explanatory text, if needed.
    * Add code samples or visuals that help complete this step.
    * Describe the result.

   Use this language:

    * _First, do x. Now, do y. Now that you have done y, do z._
    * _You must always do x before you do y because…_ followed by minimal explanation of actions in the most basic
      language possible. Link out to more detailed explanation.
    * _The output looks like this…_
    * _Notice that…/remember that…_

5. Lastly, reiterate the outcome the users will have achieved upon completion of this tutorial:
    * _You have [built, created, established, etc.]... that ..._

## How-To guides

A How-To guide instructs users on how to complete a common task. A How-To guide is always result-oriented. The goal is
to help users achieve their task in the shortest time and with the best results.

A How-To guide is like a recipe. Users don't need to know _all_ the ways in which they can complete a task. They only
need to know one optimal solution. Similarly, when users are focused on the current task, they don't want to be
distracted by the technical aspects of how things work, descriptions of concepts, or the basics of working with UI.

Avoid linking to concepts in the middle of a How-To. Every time a user clicks a link, a new tab opens. Clicking more
links opens more tabs and adds more distraction and can be overwhelming.

A How-To guide should be in this format:

* First, briefly state the goal.
* Then, list prerequisites, if any.
* Next, provide the steps to complete the task. Indicate dependencies.
* Lastly, describe the outcome and next steps (both are optional).

For the header, use the structure _active verb + noun_ (such as _Create a connection_). A troubleshooting How-To has a
specific header format:

* Include at least a partial error message. If you don't put the full error in the title, include it in the body text.
* Use fewer than 70 characters.

## Topics

Topics typically introduce concepts (such as a platform feature or capability) and answer these questions:

* What is this?
* Why is it relevant, and why would I use it?

Think of everything a user might want to know if they’ve never heard of this concept before. Notice there’s no _how_ in
these questions. Tell users _what_ it is and _why_ they need it. Reserve _how to do something_ for tutorials and how-to
guides. Link to reference material where appropriate. Don't describe more than one concept per article. Start a new
topic and link to it instead.

Topics should be in this format:

1. First, explain what this thing is.
1. Then, describe why or what users need this thing for.
1. Link to other topics and reference material when appropriate.
1. Link to How-To guides and tutorials.

For topic headers, use a noun (like _Widgets_). If a noun is ambiguous, you can add a gerund. For example, _Documenting
versions_ instead of _Versions_. Be specific:

* Either use a specific noun or a phrase that someone would search for.
* Where appropriate, use a noun followed by _workflow_. For example, _Merge request workflow_.

## Reference material

Reference material is the place to describe "what's under the hood" and provide instructions on its use. Assume that
users already understand the basic concepts involved but need to know technical details.

Unlike topics, reference material may instruct users how to do something. However, stay tightly focused on the subject
and avoid general explanations. Link to a topic that explains the concept instead. Use dry, simple language:

* State the facts: _X is an example of y. W needs to be activated using z._
* List commands, options, features, etc.: _Your options are: a,b,c, and d._
* Use warning and provide limitations: _You must have privileges to access a. Never d. Don't use b unless c._

Use tables and lists where relevant.

Like with topics, reference material headings are usually nouns.
