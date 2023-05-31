---
title: "Metadata"
---

Any [entity](../datagrok/objects.md) can be annotated with metadata, which can be used for searching. In Datagrok, there
are three kinds of metadata: tags, parameters, and properties.

:::tip

Right-click a column's header and select **Properties...** to view associated tags and other information.

:::

## Tags

A tag is a short string that can be associated with any entity. It is convenient for grouping things
(potentially of different kings) together.

To search for tagged entities, enter '#tagname' in the search box.

## Parameters

Each entity can be associated with a list of key-value pairs.

To search for entities with the specified parameters, use the `<paramName> <op> <value>` syntax, for instance:

* `imported < 10/02/2019`
* `rows > 200`
* `author = @current`

## Properties

Properties are somewhat similar to parameters, however they might contain nested properties; also, they are not directly
editable. Some properties are common for all [entities](../datagrok/objects.md), others are specific to the entity type.

Common properties:

* ID
* name
* author
* createdOn
* updatedOn
* commentedBy
* usedBy
* starredBy

## Search

Search criteria can be combined using 'and' or 'or' operands:

* `#demo or #chem`

:::tip

You can leverage metadata as search criteria in [projects](../datagrok/project.md). You can also filter tables using [smart search](../datagrok/smart-search.md) based on various fields like ID, name, rowCount, colCount, createdOn, updatedOn, author ([user](../govern/user.md) object),starredBy ([user](../govern/user.md) object), commentedBy ([user](../govern/user.md) object).


See also

* [Search](../datagrok/smart-search.md)
* [Entities](../datagrok/objects.md)
