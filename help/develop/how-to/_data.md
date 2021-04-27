<!-- TITLE: Efficiently Work with Data -->

# Data Manipulation

## Dataframe

## Bitset

## Best Practices

<!--
- use special data structures for high-performance tasks
- formatting in visualizations
- tags in dataframes, detectors
-->

### Formatting

### Tags

When annotating dataframes and columns with [metadata](../../discover/tags.md),
take a moment to consider whether you want a tag to be shown to the user. If
metadata is attached to facilitate other parts of the data workflow, or you find
it difficult to read the value stored in a tag, it makes sense to discard it
from tooltips where people can see it. To do this, simply start a tag name with
the `.` prefix.

### Detectors

Sematic type [detectors](semantic-type-detector.md) are meant to be lightweight,
simple, and efficient. Keep this fact in mind when designing these functions. If
columns you are interested in generally contain a lot of unique values, use a
special sampling method to run your checks on a random subset of column
categories. Empty values deserve special attention: make sure you don't match
them with your semantic type as it will be confusing to find it assigned to a
column consisting solely of nulls.
