---

title: Contribute docs
pagination_prev: null
---

## Markdown and GitActions

We write our technical documentation using [GitHub-flavored Markdown](https://github.github.com/gfm/) syntax. The source code and READMEs are in our [public repository on GitHub](https://github.com/datagrok-ai/public/tree/master). For our wiki, we use [Docusaurus](https://docusaurus.io/), which converts [GitHub Flavored Markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) into HTML. For basic information about writing in Docusaurus Markdown, see our [Markdown writing guide](markdown.md). To set up the environment and learn about our GitHub policy, see the [GitHub workflow page](github-workflow.mdx).

## How we think about documentation

Good documentation is clear, simple, and complete. To help you write clearly and concisely, we've created a [writing style guide](writing-style.md). We also follow these principles:

* Every page is page one. This means users can fully understand a topic and achieve their objectives by reading a single page. To follow this principle, you should provide comprehensive information, clear explanations, and relevant examples within each page. To avoid creating lengthy, overwhelming documents, we leverage [plain language](writing-style.md) and [Docusaurus Markdown features](markdown.md#content). 
* How-To pages are reserved for complex multi-step processes that require long, detailed explanations. For simpler tasks, such as launching a tool from the **Top Menu**, we use the [details feature](markdown.md/#content) or write them out in a single sentence. This ensures that users can quickly grasp and execute the task without needing to switch between different pages.
* Rather than creating separate tutorial documents, we incorporate tutorials directly into the Datagrok platform by [coding them](../how-to/write-tutorials.md). We then link to these tutorials within our documents as needed. For example, [here's a cheminformatics tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Cheminformatics/VirtualScreening). The same approach applies to code samples like scripts and functions.
* Reference pages primarily serve as lookup resources and often present information in table format. [Supported file formats](../../access/files/supported-formats.md) is a good example.  

## Before you commit to master

Before committing to master, always run these checks:

* Check for spelling or grammatical errors. Use tools like [Grammarly](https://app.grammarly.com/) or generative AI to spot and/or fix grammar, punctuation, and style errors. 
* Check that you've included the [front matter](markdown.md#metadata-front-matter).
* Check that you don't have an `H1` within the body of your doc.
* Check links and ensure they are not returning `404`.
* Check for broken images.
* Check for broken (poorly rendered) code snippets.

Once you've finished checking the basics, do the following:

* Look for opportunities to cross-link to other pages in docs.
* Focus on the doc itself and the user experience it creates:
  * Identify your target audience and ensure they can understand the content clearly. Remove any unnecessary information or details that are common knowledge and require no explanation.
  * Assess the flow of the document when a user starts reading from the beginning. Ensure that the document's structure makes it easy to understand the purpose and topic of the page.
  * Consider users who get to a specific section within the page from another page or a search engine. Evaluate whether there is enough context provided for them to understand both the section contents, and how it relates to the rest of the document.

Lastly, make sure you adhere to our [commit message policy](../dev-process/git-policy.mdx#commit-message-policy).

## Resources

If you'd like to learn more about technical writing, check out these resources:

* [Google Developer Documentation Style Guide](https://developers.google.com/style) and [Google Technical Writing Courses](https://developers.google.com/tech-writing). These resources may use a style different from ours, but still are great way to get started.
* [Every Page is Page One](https://everypageispageone.com/the-book/). A helpful overview of how to think of what goes into a page.
* [I'd rather be writing](https://idratherbewriting.com/). Guides and thoughts on tech writing process and content.
* [Write the Docs](https://www.writethedocs.org/guide/writing/beginners-guide-to-docs/). A global community that have collected a set of talk videos, articles, links, and resources.
* [Plain language resources](https://www.plainlanguage.gov/resources/) from the United States government. 