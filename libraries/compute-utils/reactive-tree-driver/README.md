# Reactive Tree Driver

Reactive tree driver, propagates data through dynamically created and
mutated funcalls trees.

## Tree structure

Tree is composed either of nested trees or FuncCallNode leaf nodes. For example:

```
[FuncCallNode, [FuncCallNode, [FuncCallNode, FuncCallNode]], FuncCallNode]
```

Here an array represents a tree level, nested arrays are nested
trees. Each tree level has one of the folowing predefined types:

- `static`, an immutable sequence (nested items could be mutable).
- `parallel`, a list of independent subtrees or FuncCallNodes, which
  could be added or removed in the runtime.
- `sequential`, a mutable sequence, where each item could use data
  from previous items.

Tree configuration is used to describe which nested items could be
added in parallel/sequential levels.

## Links

Data is propagated via links. Links are using queries to match against
the current tree state FuncCallNode inputs and outputs. Handler
defined with the link are used to transform inputs and set outputs.
Multiple link instances are allowed, using base path expand operation.

There are several roles of links:

- `data` - propagates values between FuncCalls inputs/outputs.
- `validator` - propagates validations result for inputs/outputs.
- `meta` - propagates metadata for visual ui element and visual hooks.

## Link current limitations

The limitations are not checked rn, and will cause a misbehaviour in
the runtime.

- no cycles in data links
- forward or self only based on dfs FuncCallNode ordering
- no transitive depencies in FuncCallNode inputs/outputs for data
  links (fine for validator and meta)

Also implicit backward changes are not supported. For example, the
following mutataion happend in sequential level:

```
[A, B] -> [A, B, C]
```

`C` was added at the end. If there is a link `L` that is matching only
on the new tree state and has outputs both in B, C, it will misbehave,
since adding a new step is changing the previous one via running a
newly created link. If a link `L` was matching the state before adding
`C`, it would work properly.
