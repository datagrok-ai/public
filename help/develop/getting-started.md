<!-- TITLE: Getting Started with Building on Datagrok -->
<!-- SUBTITLE: -->

# Getting Started with Building on Datagrok

Welcome to [Datagrok](https://www.datagrok.ai), the next generation data platform.

This is a short guide to jump-start with extending Datagrok and building applications on top of it.

In case of questions during the process, don't hesitate to post a question at our
[Community Forum](https://community.datagrok.ai/)! We usually respond within less than a day.

*Table of Contents:*

1. [F.A.Q. at the beginning](#faq-frequently-asked-questions)
2. [Video Walkthrough](#datarok-video-walkthrough)
3. Datagrok Tutorials ([read here first](#datagrok-tutorials-and-exercises), [start directly](https://public.datagrok.ai/))
4. Datagrok Exercises ([read here first](#datagrok-tutorials-and-exercises), [start directly](exercises.md))
5. Extending and customizing Datagrok ([read here first](#extending-and-customizing-datagrok), [start directly](https://datagrok.ai/help/develop/extending-and-customizing))
6. Building Datagrok Applications ([read here first](#datagrok-application-development), [start directly](develop/how-to/build-an-app.md))

## F.A.Q.: Frequently Asked Questions

Though these topics are outlined in a few places further, we decided to bring them at the very start.

1. _How to set up an IDE (VS Code or WebStorm) for debugging on Datagrok? Debugging doesn't work for me now._  
   The detailed guides are here:
   [VS Code](develop/develop.md#one-click-debugging-with-visual-studio-code), 
   [WebStorm](develop/develop.md#one-click-debugging-with-jetbrains-ides).  
   Follow them in detail and the debugging should work.
   Post a question with details [here](https://community.datagrok.ai/), if, afterwards, it doesn't.

2. _I've installed `npm` and `datagrok-tools` packages, but `grok` tools aren't available on my paths._  
   It's likely because the `-g` key was forgotten when installing `datagrok tools`. This key makes the
   installed package with the tool globally available for all npm packages. Re-install `datagrok-tools`
   as this: `npm install -g datagrok-tools`.

## Datagrok Video Walkthrough

Out of our diverse set of [recorded user meetings and tutorials](https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg)
we've curated those letting you start fast with building on Datagrok.
Beside getting an overall knowledge about the platform, you will directly start with the JavaScript
API and platform extensions. There's about 5 hours of watching to begin with.

It's up to you if you watch everything first and then proceed to exercises, or prefer to interleave these.

* Start with a [Coffee Company](https://www.youtube.com/watch?v=tVwpRB8fikQ) feature demo

* Proceed with all videos from [the Platform Overview](video-contents.md/#getting-started)

* Proceed with all videos [the JS API Overview](#java-script-api)

* Learn about some tooling:  
[Datagrok Tools](https://www.youtube.com/watch?v=zVVmlRorpjg&t=258s)  
[VS Code Integration](https://www.youtube.com/watch?v=zVVmlRorpjg&t=870s)  

* Learn about data access, at least skim through [some parts](video-contents.md/#data-access)

* Learn about visualizations, try to go through [all types of viewers](video-contents.md/#visualization)

* Learn about predictive modeling, at least skim through [some parts](video-contents.md/#predictive-modeling)

* Learn how our users build on top of the platform:  
[ChaRPy: converting charts to R and Python scripts](https://www.youtube.com/watch?v=seAgx5TbrzI&t=162s)
[MedChem Browser (Andrey Santrosyan & Dmitrii Petrov, GNF)](https://www.youtube.com/watch?v=seAgx5TbrzI&t=970s)

Skip all the Q&A parts at the first round. You may return to them later.

We provide a diverse set of [recorded overviews](video-contents.md) and demos for the Datagrok platform
on our [Youtube channel](https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg).
Please subscribe to this channel to receive the latest updates!

Check the full [Video Table of Contents](video-contents.md).

## Datagrok Tutorials and Exercises

It's time to take some hands-on experience with building Datagrok.

Before going to programming exercises, try the guided tutorials. They will take you step by step through
some of the platform aspects, annotating the process in the UI. Go to the `Help` button on the Datagrok's
left sidebar and open `Tutorials`. Try at least 2-3 of them.

We've prepared exercises to accommodate you with most imporant parts of Datagrok development.
The exercises tells a story of extending Datagrok with some simple bioinformatics functionality
related to nucleotide sequences.

Find the exercises here:

[https://datagrok.ai/help/develop/exercises](exercises.md)

We recommend the following:

* work on the exercises in their predefined order, as they are logically connected
* create a standalone private repository on Github under your account; this will simplify follow-ups:
   when you have questions, you may let Datagrok team members in this repo to help guide you
* rather do them 1-2 every day during the week instead of making all at once, this will help
  new concepts settle in mind and some insights pop up

## Extending and customizing Datagrok

After some hands-on experience, it's time to learn about rich Datagrok capabilities in customization.
You'd be impressed with how extendable and configurable Datagrok platform is.

Read this intro article to learn about extending and customizing Datagrok:

[https://datagrok.ai/help/develop/extending-and-customizing](extending-and-customizing.md)

## Datagrok Application Development

It's time to learn about one of the major aspects of Datagrok — building applications.
It's very likely what you'd be helping your customers with on Datagrok.

Here's a standalone how-to guide on application development:

[https://datagrok.ai/help/develop/how-to/build-an-app](develop/how-to/build-an-app.md)

Take your time to read it through and align the new knowledge with the examples linked from the guide.
Pay attention to spin-off links in this guide leading to some other guides, like these about
custom views or database connections.

## F.A.Q.: Advanced Topics

Some topics are of particular interest while learning Datagrok.

1. _Others don't see my published package on the selected Datagrok instance. How can I fix this?_  
   It is probably because you've published it under a `Debug` mode, using `grok publish`.
   In this mode, the newly published package is only seen by you. Publish it with
   `grok publish --release` to update the actual package shared to others.
   