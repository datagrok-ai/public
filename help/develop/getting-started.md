<!-- TITLE: &#8203;Getting started-->
<!-- SUBTITLE: -->

# Getting started

Welcome to Datagrok, the next-generation data platform! This document provides a step-by-step instruction to help you
start using Datagrok as an engineer.

## 1. Datagrok introduction

Our YouTube channel has a [Datagrok demo] that shows the key capabilities of the platform. Check it out for a general
understanding of Datagrok.

## 2. JavaScript API examples

Datagrok has a [a few JavaScript samples] that use the [powerful JavaScript API]. You will use the API to build
packages, applications, or plugins for the platform.

Here are a few examples of the said packages and plugins:

* [Digital Signal Processing]
* [Biosignals]
* [Cheminformatics]
* [Cell renderers]
* [Viewers]
* [Filters]

## 3. Tooling for local development

Datagrok comes with the `datagrok-tools` npm package to help you scaffold a package, an application, or a project on
your local environment. You will also use `datagrok-tools` to build and publish your code. For more information, check
out the [datagrok-tools README].

```shell
npm install --global datagrok-tools
```

<!-- PS: You might want to start by creating a [package]. -->

<!--
## Tutorial // a short tutorial that explains how to build something simple with Datagrok JavaScript API

Check out our practical tutorial that will teach you the basics of using the JavaScript API. In the tutorial, you will
generate a package, add simple code that uses the Datagrok JavaScript API, and publish this package to Datagrok's
development server.
-->

## 4. JavaScript exercises

To gain more experience with the Datagrok JavaScript API, you might want to work on a series of [JavaScript exercises].
They will introduce a more advanced usage of the API and help better understand the development process with Datagrok.

## 5. Useful links

Check out the following spaces and documents for more information about Datagrok:

1. [Community forum] and [meetings]
   To join a meeting, go to the [Community Forum] and ask for a Zoom link.
2. [Release notes]
3. [Architecture]
4. [Extending and customizing Datagrok]
5. [Building an application]
6. [Building a UI]
7. [Performance]

## 6. Videos

We've curated a few video records that introduce the platform, show the usage of the JavaScript API, and talk about
extensions.

1. [Demo]
2. [Platform overview]
3. [JS API overview]
4. [datagrok-tools overview (part 1)], [datagrok-tools overview (part 2)]
5. [VS Code Integration]
6. [Data access]
7. [Visualization and viewers]

[community forum]: https://community.datagrok.ai/

[meetings]: https://www.youtube.com/watch?v=p7_qOU_IzLM

[Release notes]: https://datagrok.ai/help/develop/release-history

[architecture]: admin/architecture.md

[Extending and customizing Datagrok]: packages/extensions.md

[Building an application]: how-to/build-an-app.md

[Building a UI]: advanced/ui.md

[Performance]: advanced/performance.md

[Demo]: https://www.youtube.com/watch?v=tVwpRB8fikQ

[Platform overview]: ../video-contents.md#getting-started

[JS API overview]: ../video-contents.md#javascript-api

[datagrok-tools overview (part 1)]: https://www.youtube.com/watch?v=zVVmlRorpjg&t=258s

[datagrok-tools overview (part 2)]: https://www.youtube.com/watch?v=0QxzllnBreI&t=4657s

[//]: # ([Building a UI]: ./ui.md)

[VS Code Integration]: https://www.youtube.com/watch?v=zVVmlRorpjg&t=870s

[Data access]: ../video-contents.md#data-access

[Visualization and viewers]: ../video-contents.md#visualizations

[Datagrok demo]: https://www.youtube.com/watch?v=tVwpRB8fikQ

[a few JavaScript samples]: https://public.datagrok.ai/js

[powerful JavaScript API]: https://datagrok.ai/js-api/

[Digital Signal Processing]: https://github.com/datagrok-ai/public/tree/master/packages/DSP

[Biosignals]: https://github.com/datagrok-ai/public/tree/master/packages/BioSignals

[Cheminformatics]: https://github.com/datagrok-ai/public/tree/master/packages/Chem

[//]: # ([Natural Language Processing]: https://github.com/datagrok-ai/public/tree/master/packages/NLP)

[Cell renderers]:https://github.com/datagrok-ai/public/blob/master/libraries/chem-meta/src/rdkit-api.ts

[viewers]: https://github.com/datagrok-ai/public/tree/master/packages/Viewers

[filters]: https://github.com/datagrok-ai/public/blob/master/packages/Widgets/src/filters/radio-button-filter.ts

[datagrok-tools README]: https://github.com/datagrok-ai/public/tree/master/tools#datagrok-tools

[JavaScript exercises]: exercises/exercises.md
