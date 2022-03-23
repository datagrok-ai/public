<!-- TITLE: Metadata -->
<!-- SUBTITLE: -->

# Metadata

Any [entity](../overview/objects.md) can be annotated with metadata, which can be used for searching. In Datagrok, there
are three kinds of metadata: tags, parameters, and properties.

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
editable. Some properties are common for all [entities](../overview/objects.md), others are specific to the entity type.

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

See also

* [Search](../overview/smart-search.md)
* [Entities](../overview/objects.md)
