<!-- TITLE: &#8203;Getting Started-->
<!-- SUBTITLE: -->

# Getting Started with Building on Datagrok

Welcome to [Datagrok](https://www.datagrok.ai), the next generation data platform!

This is a guide to help you jump-start extending Datagrok via plugins and building applications on top of it.

During this Getting Started, don't hesitate to post questions and suggestions at our
[Community Forum](https://community.datagrok.ai/). We usually respond within less than a day.

*Table of Contents:*

1. [Daily Sources of Information](#daily-sources-of-information)
2. [Datagrok Video Walkthrough](#datarok-video-walkthrough)
3. [Quick F.A.Q., before you start developing](#faq-frequently-asked-questions)
4. Datagrok Tutorials ([read here first](#datagrok-tutorials-and-exercises), [start directly](https://public.datagrok.ai/))
5. Datagrok Exercises ([read here first](#datagrok-tutorials-and-exercises), [start directly](exercises.md))
6. Datagrok Architecture ([read here first](#datagrok-architecture), [start directly](develop/admin/architecture.md))
7. Extending and customizing Datagrok ([read here first](#extending-and-customizing-datagrok), [start directly](https://datagrok.ai/help/develop/extending-and-customizing))
8. Building Datagrok Applications ([read here first](#datagrok-application-development), [start directly](develop/how-to/build-an-app.md))

## Daily Sources of Information

Recall the following locations during your daily work with Datagrok:

1. We post regular updates to our [Community Forum](https://community.datagrok.ai/) on the following topics:

   * [JavaScript API Updates](https://community.datagrok.ai/t/javascript-api-updates/526/9)
   
   * [Visualization-related Updates](https://community.datagrok.ai/t/visualization-related-updates/521/12)
   
   * [Cheminformatics Updates](https://community.datagrok.ai/t/cheminformatics-updates/457/9)
   
   We recommend you to subscribe to at least these three topics.
   
2. We aggregate Datagrok Release Notes [here](https://datagrok.ai/help/develop/release-history).
   For every build, the notes are auto-generated. However, for the major releases we [curate and shape them up](https://github.com/datagrok-ai/public/blob/master/help/develop/release-history.md#2021-04-14-dev-build-08936)
   by hand with the lists of highlights, major items, improvements, bug fixes and cross-linking to relevant Community Forum posts,
   examples and sources. Therefore, Notes become a powerful source of full knowledge about what's coming in the recent Datagrok update.

3. Highly diversed set of [JS API Samples](https://public.datagrok.ai/js) right at hand. Use the built-in
   [JS Fiddler](https://public.datagrok.ai/js) inside Datagrok to fast-prototype interfaces and functions
   which will later become parts of your packages. The sources for these API Samples are versioned
   [here](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples).
   
4. [Datagrok Public repository](https://github.com/datagrok-ai/public). We are passionate about making
   more parts of our work open-source, as well as on involving other developers into expanding the platform.
   The results of this you may find on our [public github](https://github.com/datagrok-ai/public).  
   
   Not only this is a great source of really many API use examples, but also a repository of fine-crafted
   patterns, such as a [Chembl data repository browser](), and plugins, such as
   [custom cell renderers](https://github.com/datagrok-ai/public/blob/master/packages/Chem/src/rdkit_cell_renderer.js),
   [custom viewers](https://github.com/datagrok-ai/public/tree/master/packages/Viewers),
   [custom filters](https://github.com/datagrok-ai/public/blob/master/packages/Widgets/src/filters/radio_button_filter.js), and so on. Also find there our packages for domains such as
   [Digital Signal Processing](https://github.com/datagrok-ai/public/tree/master/packages/DSP), [Biosignals](https://github.com/datagrok-ai/public/tree/master/packages/BioSignals), [Cheminformatics](),
   or [Natural Language Processing](https://github.com/datagrok-ai/public/tree/master/packages/NLP).

5. We host regular Datagrok User Meetings. Ask Andrew or Dan for a Zoom link
   (for example, by a private message at [Community Forum](https://community.datagrok.ai/)).

## Datagrok Video Walkthrough

Let's get it started!

From our detailed range of [recorded user meetings and tutorials](https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg)
we've curated these which let you start fast with building on Datagrok. Besides getting an overall idea about the platform,
which will help you a lot all the way further in your work, in these videos you'll be directly introduced to our
JavaScript API and platform extensions.

There's about 5 hours of watching to begin with. It's up to your taste if you watch everything first and then
proceed to [exercises](#datagrok-tutorials-and-exercises), or prefer combining them on the way.

1. Start with a [Coffee Company](https://www.youtube.com/watch?v=tVwpRB8fikQ) Feature Demo to get an overall idea

2. Going further, proceed with all videos from [the Platform Overview](video-contents.md/#getting-started)

3. Using Datagrok JS API would become your regular routine. Proceed with all videos [the JS API Overview](video-contents.md/#java-script-api)

4. Learn about some daily Datagrok tooling:  

   * [Datagrok Tools](https://www.youtube.com/watch?v=zVVmlRorpjg&t=258s)
   
   * [VS Code Integration](https://www.youtube.com/watch?v=zVVmlRorpjg&t=870s)

5. Datagrok is well tailored for data integration. Learn about data access, at least skim through [some parts](video-contents.md/#data-access)

6. Visualization is a crucial part of Datagrok. Learn about visualizations, try to go
   through [all types of viewers](video-contents.md/#visualization)

7. Learn about Datagrok's predictive modeling, at least skim through [some parts](video-contents.md/#predictive-modeling)

8. Get a taste of building on top of the platform:  

   * [ChaRPy: converting charts to R and Python scripts](https://www.youtube.com/watch?v=seAgx5TbrzI&t=162s)
   
   * [MedChem Browser (Andrey Santrosyan & Dmitrii Petrov, GNF)](https://www.youtube.com/watch?v=seAgx5TbrzI&t=970s)

Skip all the Q&A parts at the first round. You may want to return to them later.

All the [recorded overviews](video-contents.md) and demos for the Datagrok platform  are on our [Youtube channel](https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg). Please subscribe to it and receive the latest updates!
[Here is](video-contents.md) the current videos' [full Table of Contents](video-contents.md).

## F.A.Q.: Frequently Asked Questions

We've learned from experience it's very likely you'd have these questions at the beginning. So, we curated them for you in advance.

*Q:* How to set up an IDE (VS Code or WebStorm) for debugging on Datagrok? Debugging doesn't work for me now.  

*A:* The detailed guides are here: [VS Code](develop/develop.md#one-click-debugging-with-visual-studio-code),
[WebStorm](develop/develop.md#one-click-debugging-with-jetbrains-ides).   Post a question with details [here](https://community.datagrok.ai/), if debugging won't work after taking the above steps.

*Q:* I've installed `npm` and `datagrok-tools` packages, but `grok` tools aren't available on my paths.  

*A:* It's likely because the `-g` key was forgotten when installing `datagrok tools`. This key makes the
     installed package with the tool globally available for all npm packages and locations. Re-install `datagrok-tools`
     as this: `npm uninstall datagrok-tools` and `npm install -g datagrok-tools`.

## Datagrok Tutorials and Exercises

It's time to get into hands-on experience with building Datagrok.

Before going to programming exercises, try the guided tutorials. They will take you step by step through
some of the platform aspects, annotating the process in the Datagrok's UI. Go to the `Help` button
on the Datagrok's left sidebar and open `Tutorials`. Try at least 2-3 of them.

We've prepared exercises to accommodate you with the most important parts of Datagrok development.
The exercises tells a story of extending Datagrok with some simple bioinformatics functionality
related to nucleotide sequences. The exposed functionality itself is quite schematic to not steal
too much of your time for each exercise. However, the use of Datagrok API is all based on real-world
cases.

Find the exercises here:

[https://datagrok.ai/help/develop/exercises](exercises.md)

We recommend the following:

* pay attention to the first part with setting up the environment
* work on the exercises in their predefined order, as they are logically connected
* create a standalone private repository on Github under your account; this will simplify follow-ups —
  when you have questions, you may let Datagrok team members in this repo to help guide you
* though it's a lot about personal preferences, we recommend doing 1-2 exercises every day during the week
  instead of making them all at once; this will help new concepts settle in mind and some insights pop up
  
## Datagrok Architecture
 
We are highly interested in helping developers building on top of Datagrok deliver the best experience
to their end-users. That's why we highly encourage you to go through these two articles on architecture,
either before or after doing the Exercises:

[https://datagrok.ai/help/develop/admin/architecture](develop/admin/architecture.md)

[https://datagrok.ai/help/develop/performance](develop/performance.md)

There are architectural pillars regarding accessing data and leveraging the power of our
in-memory columnar data store. This understanding is crucial for making your applications
and plugins low-latency and smooth in operation.

## Extending and customizing Datagrok

After some hands-on experience, it's time to learn about rich Datagrok capabilities in customization.
You'd be impressed with how extendable and configurable Datagrok platform is!

Read this intro article to learn about extending and customizing Datagrok:

[https://datagrok.ai/help/develop/extending-and-customizing](extending-and-customizing.md)

## Datagrok Application Development

It's time to learn about one of the major aspects of Datagrok — building applications.
It's very likely what you'd be helping your customers with building their solutions on Datagrok.

Here's a standalone how-to guide on Datagrok application development:

[https://datagrok.ai/help/develop/how-to/build-an-app](develop/how-to/build-an-app.md)

In a logical sequence it hilights the aspects important to this subject. Take your time to read it through
and align the new knowledge with the examples linked from the guide. Pay attention to spin-off links in this guide
leading to some other guides, like these about [custom views](develop/how-to/custom-views.md) or
[database connections](develop/how-to/access-data.md).

## F.A.Q.: Advanced Topics

Some topics are of particular interest while learning Datagrok.

*Q.* Others don't see my published package on the selected Datagrok instance. How can I fix this?  

*A.* It is probably because you've published it under a `Debug` mode, using `grok publish`.
     In this mode, the newly published package is only seen by you. Publish it with
     `grok publish --release` to update the actual package shared to others.
   