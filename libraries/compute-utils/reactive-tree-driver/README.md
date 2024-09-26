# Reactive Tree Driver

Reactive tree driver, propagates data through dynamically created and
mutated funcalls trees.

## Tree structure

Tree is composed either of nested trees or FuncCallNode leaf
nodes. For example:

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

A tree configuration is used to describe which nested items could be
added on parallel/sequential levels.

## Links

Data is propagated via links. Links are using queries to match against
the current tree state FuncCallNode inputs and outputs. A handler
funciton defined or referenced with links is used to transform inputs
and set outputs.  Multiple link instances are allowed, using base path
expand operation.

There are several roles of links:

- `data` - propagates values between FuncCalls inputs/outputs.
- `validator` - propagates validations result for inputs/outputs.
- `meta` - propagates metadata for visual ui element and visual hooks.

### Links current limitations

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

### Links tree mutations rerunning

When a tree is mutated links are recreated and selectivly rerun. For
example, node `I` was added, nodes here could be nested subtrees:

```
[A, [B, C, D], E] -> [A, [B, C, I, D], E]
^	^                ^	 ^
|   |                |   |
X   Y                X   Y
```

- Data links that has outputs either in `I` or after on the same level
and inputs before `I` will be rerun. If there is an output before `I`,
it will be discarded.
- Data links the has inputs in `I` and outputs outside of `I` will be
  run.
- Meta links that have either input or output in `Y` will be rerun.
- All validators will be rerun.
- All orphaned metadata will be discarded.

Note that links run order is based on the DFS order of the first
encountered input. After data link is triggered, all triggered
dependent links are awaited, then the next one is executed.


For there removing case:
```
[A, [B, C, I, D], E] ->  [A, [B, C, D], E]
^	^                    ^   ^
|   |                    |   |
X   Y                    X   Y
```

- Data links the had outputs after deleted element `I` (index 2 from
0) on the same level and inputs before will be rerun. If there is an
output before, it will be discarded.
- Data links the has inputs in `Y` and outputs outside of `Y` will be
  run.
- Meta links that have either input or output in `Y` will be rerun.
- All validators will be rerurn.
- All orphaned metadata will be discarded.

Similar rules are applied to the moving case.
