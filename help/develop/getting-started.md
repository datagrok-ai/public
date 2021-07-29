<!-- TITLE: &#8203;Getting Started-->
<!-- SUBTITLE: -->

# Getting started with building on Datagrok

Welcome to Datagrok, the next-generation data platform!

This guide will help you jump-start using Datagrok as an engineer. In short, what you can do with Datagrok is extend it
with plugins, build your applications on top of Datagrok, perform data science tasks, run scripts, retrieve and explore 
data, and many other things.

If you have any questions or suggestions after reading this Getting Started guide, go to our [Community Forum] &mdash; 
we typically respond within 24 hours.

## Let's Start

If you already know something about Datagrok (probably thanks to our [YouTube channel]), you can take a fast track:

* Have a look at the [20-minute introduction to a Datagrok application] for a refresher, and then proceed to the
  [JavaScript programming exercises](exercises.md) that explain how to work with Datagrok using the JavaScript API.
* Check out the [Daily Sources of Information](#daily-sources-of-information), where you can post questions and get
  information about the platform updates. Return to Daily Sources of Information from time to time to study Datagrok in
  depth.

In case you want a clean start, review the following sections one by one, in order:

1. [Daily Sources of Information](#daily-sources-of-information)
2. [Datagrok Video Walkthrough](#datagrok-video-walkthrough)
3. [Frequently Asked Questions](#frequently-asked-questions)
4. [Datagrok Tutorials](#datagrok-tutorials-and-exercises) or [start using the public version of Datagrok]
5. [Datagrok Exercises](#datagrok-tutorials-and-exercises) or [start directly](exercises.md)
6. [Datagrok Architecture](#datagrok-architecture) or [start directly](admin/architecture.md)
7. [Extending and Customizing Datagrok](#extending-and-customizing-datagrok) or [start extending and customizing Datagrok]
8. [Building Datagrok Applications](#datagrok-application-development) or [start directly](how-to/build-an-app.md)

## Daily sources of information

Check out the following spaces for more information about Datagrok:

1. [Community Forum]. This is where we post updates about the Datagrok platform and ecosystem. We recommend that you
   subscribe at least to the following three topics on the forum:
   * [JavaScript API Updates]
   * [Visualization-related Updates]
   * [Cheminformatics Updates]

2. [Datagrok Release History] &mdash; this page aggregates the release notes. For every new Datagrok build, the release notes
   are automatically generated. However, for the major releases, we curate our release notes, which you can see on
   [our GitHub page]. There, we detail the highlights, major items, improvements, and bug fixes in the releases. We also
   provide links to the relevant [Community Forum] topics, examples, and sources.

3. [JS API Samples]. We add diverse examples of our JavaScript API to the Datagrok public website. You can use the
   built-in JS Fiddler inside Datagrok to rapidly prototype interfaces and functions that will later become parts of
   your packages.

   The sources of these API samples are stored in the [Datagrok ApiSamples package on GitHub].

4. [Datagrok public repository]. We're passionate about making more parts of our work open-source, and we strive to
   involve other developers into expanding the Datagrok platform.

   The Datagrok public repository is not only a great source of API examples, but also a repository of fine-crafted
   patterns, such as a [ChEMBL data repository browser], and plugins, such as [custom cell renderers], [custom viewers],
   [custom filters], and others.

   Finally, you can also find our packages for domains such as [Digital Signal Processing], [Biosignals], [Cheminformatics],
   and [Natural Language Processing].

5. Datagrok User Meetings. We host regular meetings with active Datagrok users. If you want to join such meetings, go
   to the [Community Forum] and ask Andrew or Dan for a Zoom link.

## Datagrok video walkthrough

[Our YouTube channel] contains a lot of video tutorials and meetings with Datagrok users, so our channel is the best place
to start learning Datagrok. We've curated a few records that will give you an overall idea of the platform and introduce 
to our JavaScript API and various platform extensions.

1. [Coffee Company demo] &mdash; this video will give you an overall idea of Datagrok.

2. [Platform overview]

3. [JS API Overview]

4. Datagrok tooling:

   * [Datagrok Tools]

   * [VS Code Integration]

5. [Data access] &mdash; Datagrok is greatly tailored for data integration.

6. [Visualization and viewers]. Visualization is a key part of Datagrok.

7. [Predictive modeling]

8. Get a taste of building on top of the platform:

   * [ChaRPy: Converting charts to R and Python scripts]

   * [MedChem Browser (Andrey Santrosyan & Dmitrii Petrov, GNF)]

A little piece of advice: When viewing the videos, skip the Q&A parts. You might want to review those parts later.

All the overviews and demos for the Datagrok platform are located on [our YouTube channel] &mdash; subscribe to this
channel to get updates on the latest Datagrok features.

## Datagrok tutorials and exercises

It's time to get your hands dirty with building on Datagrok!

First, try out a few tutorials to learn the key Datagrok aspects. To open the tutorials, click **Help** on the Datagrok's
left sidebar and then select **Tutorials**.

Besides tutorial, we've prepared [some exercises] to teach you the basics of the Datagrok API. Through the exercises, you will
extend Datagrok with a bioinformatics functionality related to nucleotide sequences (such as DNA). These exercises,
although vastly simplified, are based on real-world cases.

Here are a few recommendations when you start the exercises:

* Make sure that you configure the environment to run Datagrok
* Complete exercises one by one, in order
* Push the changes to your private repository on GitHub, which will greatly simplify our work when you reach to clarify
  some aspects of using Datagrok API

Last but not least, you might want to complete one or two exercises per day, not rush through all of them in one day.
Additional time will help you understand the new concepts.

## Datagrok architecture

The key to efficient development on Datagrok is understanding the Datagrok architecture and performance optimizations that
you can make to bring the best results to the end user.

We encourage you to go through these two articles:

* [Datagrok Architecture]
* [Datagrok Performance]

Properly accessing data and leveraging the power of our in-memory data store is crucial to make your applications
and plugins low-latency and smooth in operation.

## Extending and customizing datagrok

After some hands-on experience, you might want to learn about rich customization capabilities of Datagrok. Check out
[the Extending and Customizing Datagrok article] for more information.

## Datagrok application development

[Build a Datagrok App] is the place to start developing your first application on our platform.

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
[start using the public version of Datagrok]: https://public.datagrok.ai/
[start extending and customizing Datagrok]: https://datagrok.ai/help/develop/extending-and-customizing
[JavaScript API Updates]: https://community.datagrok.ai/t/javascript-api-updates/526/9
[Visualization-related Updates]: https://community.datagrok.ai/t/visualization-related-updates/521/12
[Cheminformatics Updates]: https://community.datagrok.ai/t/cheminformatics-updates/457/9
[Datagrok Release History]: https://datagrok.ai/help/develop/release-history
[our GitHub page]: https://github.com/datagrok-ai/public/blob/master/help/develop/release-history.md#dev-build-08936
[JS API Samples]: https://public.datagrok.ai/js
[Datagrok ApiSamples package on GitHub]: https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples
[Datagrok public repository]: https://github.com/datagrok-ai/public
[ChEMBL data repository browser]: https://www.ebi.ac.uk/chembl/
[custom cell renderers]: https://github.com/datagrok-ai/public/blob/master/packages/Chem/src/rdkit_cell_renderer.js
[custom viewers]: https://github.com/datagrok-ai/public/tree/master/packages/Viewers
[custom filters]: https://github.com/datagrok-ai/public/blob/master/packages/Widgets/src/filters/radio_button_filter.js
[Digital Signal Processing]: https://github.com/datagrok-ai/public/tree/master/packages/DSP
[Biosignals]: https://github.com/datagrok-ai/public/tree/master/packages/BioSignals
[Cheminformatics]: https://datagrok.ai/help/domains/chem/cheminformatics
[Natural Language Processing]: https://github.com/datagrok-ai/public/tree/master/packages/NLP
[Our YouTube channel]: https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg
[Coffee Company feature demo]: https://www.youtube.com/watch?v=tVwpRB8fikQ
[Datagrok Tools]: https://www.youtube.com/watch?v=zVVmlRorpjg&t=258s
[VS Code Integration]: https://www.youtube.com/watch?v=zVVmlRorpjg&t=870s
[Platform overview]: ../video-contents.md/#getting-started
[JS API Overview]: ../video-contents.md/#java-script-api
[Data access]: ../video-contents.md/#data-access
[Predictive modeling]: ../video-contents.md/#predictive-modeling
[ChaRPy: Converting charts to R and Python scripts]: https://www.youtube.com/watch?v=seAgx5TbrzI&t=162s
[MedChem Browser (Andrey Santrosyan & Dmitrii Petrov, GNF)]: https://www.youtube.com/watch?v=seAgx5TbrzI&t=970s
[our Youtube channel]: https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg
[Architecture]: admin/architecture.md
[Performance]: performance.md
[some exercises]: exercises.md
[the Extending and Customizing Datagrok article]: extending-and-customizing.md
[Build a Datagrok App]: how-to/build-an-app.md
[VS Code]: develop.md#one-click-debugging-with-visual-studio-code
[WebStorm]: develop.md#one-click-debugging-with-jetbrains-ides