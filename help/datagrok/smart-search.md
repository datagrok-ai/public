<!-- TITLE: Smart search -->
<!-- SUBTITLE: -->

# Smart search

Use the free-text input that lets define complex queries. Smart search supports AND and OR operators and parenthesis, so
you can combine filters. If you type single string - search engine will treat it as filter by name. Tags filtering is
supported: #demo will show entities, tagged by #demo tag, also you can combine tags conditions using AND or OR
operators. Every entity has properties, that could be used for filtering. [See more](../datagrok/objects.md).

## Examples

Unstructured query; looks for 'biologics' in title and description:

```
Biologics
```

Having #demo tag:

```
# Demo
```

Tagged as either either #demo or #chem:

```
# Demo or #chem
```

Created in the last 7 days:

```
createdOn > -1w
```

Complex conditions:

```
(#demo and #chem) or author = "john@google.com"
starredBy = @current or author = @current
```

Created by recently joined users:

```
author.joined > -5d
```

See also:

* [Find and replace](../transform/find-and-replace.md)
