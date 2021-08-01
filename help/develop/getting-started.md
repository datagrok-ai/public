<!-- TITLE: &#8203;Getting started-->
<!-- SUBTITLE: -->

# Getting started with building on Datagrok

Welcome to [Datagrok], the next-generation data platform!

This guide will help you jump-start using Datagrok as an engineer. In short, what you can do with Datagrok is extend it
with plugins, build your applications on top of the platform, perform data science tasks, run scripts, retrieve and explore 
data, and many other things.

If you have any questions or suggestions after reading this guide, go to our [Community Forum] &mdash; 
we typically respond within 24 hours.

## Let's start

If you already know something about Datagrok (probably thanks to our [YouTube channel]), you can take a fast track:

* Have a look at the [20-minute introduction to a Datagrok application] for a refresher, and then proceed to the
  [JavaScript programming exercises](exercises.md) that explain how to work with Datagrok using the JavaScript API.
* Check out the [daily sources of information](#daily-sources-of-information), where you can post questions and get
  information about the platform updates. Return to these information sources from time to time to study Datagrok in
  depth.

If you want a clean start, review the following sections one by one, in order:

1. [Daily sources of information](#daily-sources-of-information)
2. [Datagrok video walkthrough](#datagrok-video-walkthrough)
3. [Datagrok exercises](#datagrok-exercises)
4. Further reading:
   * [Extending and customizing Datagrok]
   * [Building Datagrok applications]
   * [Datagrok UI](ui.md)
   * [Datagrok architecture](admin/architecture.md)
5. [Frequently asked questions](#frequently-asked-questions)

## Daily sources of information

Check out the following spaces for more information about Datagrok.

1. [Community Forum]. This is where we post updates about the Datagrok platform and ecosystem. We recommend that you
   subscribe at least to the following three topics on the forum:
   * [JavaScript API updates]
   * [Visualization-related updates]
   * [Cheminformatics updates]

2. [Datagrok release history]. For every new Datagrok build, the  release notes are automatically generated. However,
   for the major releases, we curate our release notes, which you can see on [our GitHub page]. There, we detail the
   highlights, major items, improvements, and bug fixes in the releases. We also provide links to the relevant
   [Community Forum] topics, examples, and sources.

3. [JS API samples]. We add diverse examples of our JavaScript API to the Datagrok public website. You can use the
   built-in JS Fiddler inside Datagrok to rapidly prototype interfaces and functions that will later become parts of
   your packages.

   The sources of these API samples are stored in the [Datagrok ApiSamples package on GitHub].

4. [Datagrok public repository]. We are keen on making more parts of our work open-source. We would like to involve
   other developers into expanding the Datagrok platform via this repository.

   The Datagrok public repository is not only a great source of API examples, but also a repository of fine-crafted
   examples and plugins, such as [custom cell renderers], [custom viewers], [custom filters], and others.

   Finally, you can find our packages for domains such as [Digital Signal Processing], [Biosignals], 
   [Cheminformatics], and [Natural Language Processing].

5. [Datagrok user meetings]. We host regular meetings with active Datagrok users. If you want to join the meetings, go
   to the [Community Forum] and ask Andrew or Dan for a Zoom link.

## Datagrok video walkthrough

[Our YouTube channel] contains a lot of video tutorials and meetings with Datagrok users. It's a great place to
start learning Datagrok. Below we've curated a few records that will give you an overall idea of the platform and
introduce  to our JavaScript API and various platform extensions.

1. [Coffee Company demo] (20 m.). An overall idea of Datagrok
2. [Platform overview] (1 h. 30 m.)
3. [JS API Overview] (1 h. 20 m.)
4. Datagrok tooling:
   * Datagrok Tools (
     [part 1](https://www.youtube.com/watch?v=zVVmlRorpjg&t=258s),
     [part 2](https://www.youtube.com/watch?v=0QxzllnBreI&t=4657s)) (20 m.)
   * [VS Code Integration] (10 m.)
5. [Data access] (1 h. 30 m.). Datagrok is greatly tailored for data integration
6. [Visualization and viewers] (1 h. 30 m.). Visualization is a key part of Datagrok
7. Get a taste of building on top of the platform:
   * [ChaRPy: Converting charts to R and Python scripts] (15 m.)
   * [MedChem Browser (Andrey Santrosyan & Dmitrii Petrov, GNF)] (20 m.)

When viewing the videos, skip the Q&A parts. You might want to review those parts later.

Subscribe to our YouTube channel to get updates on the latest Datagrok features.

## Datagrok exercises

It's time to get your hands dirty with developing on the Datagrok platform!

We've prepared [some JavaScript exercises] to teach you the basics of the Datagrok API. Through the 
exercises, you will extend Datagrok with a bioinformatics functionality related to nucleotide sequences
(such as DNA). These exercises, although vastly simplified, are based on real-world cases.

Here are a few recommendations when you start doing the exercises:

* Make sure that you configure the environment to run Datagrok
* Complete exercises one by one, in order
* Push the changes to your private repository on GitHub, which will greatly simplify our work when you reach to clarify
  some aspects of using Datagrok API

Last but not least, you might want to complete one or two exercises per day, not rush through all of them in one day.
Additional time will help you understand the new concepts. It typically takes a week to complete them.

You may also find useful our interactive tutorials to familiarize yourself with the different
functionalities of the platform: click **Help** on the Datagrok's left sidebar and
then select **Tutorials**.

## Further reading

1. Learn about extending Datagrok with plugins:
   
      * [Extending and customizing Datagrok]
      
2. Learn about building UI with Datagrok:

      * [Datagrok UI](ui.md)

3. Learn about building applications on top of Datagrok:

      * [Building Datagrok applications]

4. Understanding the Datagrok architecture and performance of the underlying entities is pivotal for
   bringing the best user experience. Here is where to learn more on that:

      * [Datagrok architecture]
      * [Datagrok performance]

## Frequently asked questions

*Question*
How do I set up VS Code or WebStorm for debugging on Datagrok? Debugging doesn't work for me now.

*Answer:*
Check out our guides on [VS Code] and [WebStorm]. If you configured VS Code or WebStorm following the guides, but
debugging still doesn't work, ask our technical gurus on the [Community Forum].

*Question:*
I've installed `npm` and `datagrok-tools` packages, but `grok` tools aren't available on my paths.

*Answer:*
It's likely that the `datagrok-tools` library wasn't installed globally. Install `datagrok-tools` globally:

```sh
# First, remove the locally installed package...
npm uninstall datagrok-tools

# and then install datagrok-tools globally using the -g flag
npm install -g datagrok-tools
```

*Question:*
Others don't see my published package on the selected Datagrok instance. How can I fix this?

*Answer:*
If you publish a package in debug mode using the standard `grok publish` command, only you will see the published
package. Run `grok publish --release` so others could see your package, too.

[Datagrok]: https://www.datagrok.ai
[Community Forum]: https://community.datagrok.ai/
[YouTube channel]: https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg
[20-minute introduction to a Datagrok application]: https://www.youtube.com/watch?v=tVwpRB8fikQ
[Build an App]: how-to/build-an-app.md
[JavaScript API updates]: https://community.datagrok.ai/t/javascript-api-updates/
[Visualization-related updates]: https://community.datagrok.ai/t/visualization-related-updates/
[Cheminformatics Updates]: https://community.datagrok.ai/t/cheminformatics-updates/
[Datagrok Release History]: https://datagrok.ai/help/develop/release-history
[our GitHub page]: https://github.com/datagrok-ai/public/blob/master/help/develop/release-history.md#dev-build-08936
[JS API Samples]: https://public.datagrok.ai/js
[Datagrok ApiSamples package on GitHub]: https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples
[Datagrok public repository]: https://github.com/datagrok-ai/public
[Datagrok user meetings]: https://www.youtube.com/watch?v=p7_qOU_IzLM
[ChEMBL data repository browser]: https://www.ebi.ac.uk/chembl/
[custom cell renderers]: https://github.com/datagrok-ai/public/blob/master/packages/Chem/src/rdkit_cell_renderer.js
[custom viewers]: https://github.com/datagrok-ai/public/tree/master/packages/Viewers
[custom filters]: https://github.com/datagrok-ai/public/blob/master/packages/Widgets/src/filters/radio_button_filter.js
[Digital Signal Processing]: https://github.com/datagrok-ai/public/tree/master/packages/DSP
[Biosignals]: https://github.com/datagrok-ai/public/tree/master/packages/BioSignals
[Cheminformatics]: https://github.com/datagrok-ai/public/tree/master/packages/Chem
[Natural Language Processing]: https://github.com/datagrok-ai/public/tree/master/packages/NLP
[Our YouTube channel]: https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg
[Coffee Company demo]: https://www.youtube.com/watch?v=tVwpRB8fikQ
[VS Code Integration]: https://www.youtube.com/watch?v=zVVmlRorpjg&t=870s
[Platform overview]: ../video-contents.md#getting-started
[JS API Overview]: ../video-contents.md#javascript-api
[Data access]: ../video-contents.md#data-access
[Visualization and viewers]: ../video-contents.md#visualizations
[Predictive modeling]: ../video-contents.md#predictive-modeling
[ChaRPy: Converting charts to R and Python scripts]: https://www.youtube.com/watch?v=seAgx5TbrzI&t=162s
[MedChem Browser (Andrey Santrosyan & Dmitrii Petrov, GNF)]: https://www.youtube.com/watch?v=seAgx5TbrzI&t=970s
[our Youtube channel]: https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg
[Datagrok architecture]: admin/architecture.md
[Datagrok performance]: performance.md
[some JavaScript exercises]: exercises.md
[Extending and customizing Datagrok]: extending-and-customizing.md
[Building Datagrok applications]: how-to/build-an-app.md
[VS Code]: develop.md#one-click-debugging-with-visual-studio-code
[WebStorm]: develop.md#one-click-debugging-with-jetbrains-ides