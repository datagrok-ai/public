# Writing a blog post

Posts are published at `https://datagrok.ai/blog/<slug>`.

## Layout

One folder per post so images live next to the text:

```text
blog/
  YYYY-MM-DD-short-name/
    index.md
    cover.png
  authors.yml
  tags.yml
```

The folder date is the publish date. The URL comes from `slug` in front matter, not the folder name.

## Front matter

```yaml
---
slug: short-name            # URL: /blog/short-name. Keep it stable; it is the canonical URL.
title: Post title           # <title> and h1
description: One sentence.  # Meta description and search snippet. 150 chars or less.
authors: [datagrok]         # Keys from authors.yml. Inline authors fail the build.
tags: [product]             # Keys from tags.yml. Inline tags fail the build.
image: ./cover.png          # Social card. 1200x630 recommended.
keywords: [a, b, c]         # Meta keywords.
draft: true                 # Optional. Excluded from production builds.
---
```

## Body

- Put `<!-- truncate -->` after the first paragraph or two. The build fails without it, because the list page would otherwise duplicate the full post.
- Link to help pages with absolute site paths: `[Viewers](/help/visualize/viewers/viewers)`. The build fails on broken links.
- Link to feeds and other non-route files with `pathname:///blog/rss.xml`.
- Images: `![Alt text](./cover.png)`. Always include alt text.

## Adding an author or tag

Edit `authors.yml` or `tags.yml`. Keep tags few; every tag is a public indexable page.

## Preview locally

```bash
cd docusaurus && npm run start
```

Then open `http://localhost:3000/blog`.
