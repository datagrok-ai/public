# Chem package developer guide

* If you are a contributor to this package, be guided by this document on the library
design and engineering choices.

* Contribute to this file if you have a knowledge about third-parties such as a
chemical library (RdKit, OpenChemLib, etc.) or peculiarities of Chrome, or any other
knowledge which isn't explicitly expressed in the code

* Contributors from Datagrok should also check the videos recorded about the design,
ask @dskatov and @StLeonidas on getting access to these videos and a walk-through

* If you are an external contributor, post any question on Chem directly in
[issues](https://github.com/datagrok-ai/public/issues) with a label `Question` and
a project `Chem`

## API design

## Implementation

### Compute caching

### Renders caching

## Web Workers

### Message passing interface

### TypeScript and Webpack

### Passing raw objects

## RdKit

### Source of the WASM library

@ptosco assembles special builds of RdKit WASM library (called MinimalLib) from his fork:

https://github.com/ptosco/rdkit/tree/master/Code/MinimalLib

The reason for it is that this fork ships customer-specific things.

### Scaffold rendering offset

In a method
[_drawMolecule]( https://github.com/datagrok-ai/public/blob/ad9bbbfc10347a1947a67762c19c96f4b1a0735f/packages/Chem/src/rendering/rdkit-cell-renderer.ts#L173)
(part of `Chem/src/rdkit-cell-renderer.ts`) an `offscreenCanvas` is used first to draw a
molecule, and only then this canvas contents are brought to the actual canvas
`onscreenCanvas`, which users actually see. There is a reason to such inefficiency.

If we've drawn through an RdKit method directly to the onscreen canvas, it
would've been seen that the molecule scaffold, also drawn by RdKit, is a bit off
the origin structure. This phenomenon is avoided with using an intermediate canvas.
This problem only appears when the molecules are rendered to the Datagrok's
canvas, but isn't _yet_ reproduced "in the wild".

Note that the _drawMolecule_ is only invoked when the molecule appears in the
Datagrok platform context first time, i.e. when it is first rendered to the
screen. After the first appearance, the render is put to the `rendersCache`,
and the corresponding molecule is put to `molCache`. Only when the cache record
for the molecule is evoked, the render will be required again. Therefore,
this double-drawing only takes place on the first appearance.

We need to resolve this part eventually, but the priority isn't the highest.

### The testbed

The folder `Chem/testbed` offers a standalone web-page for testing RdKit.
This page doesn't require Datagrok in place. This testbed is useful when certain
issues of the library are identified and a bug report needs to be produced.

To run the testbed:

1. Install `npm install http-server'

2. Run `run-server.bat` and go to one of the locations from the screen.
Typically it is `https://localhost:81`.

Testbed needs to run through the web server, as it serves WebAssembly.

If you test a new version of RdKit, just replace the locations in the
`GettingStartedJS.html` to the new ones and commit these changes along
with the RdKit library version being tested.

The test form is forked from RdKit's original testbed:
[link](https://github.com/ptosco/rdkit/tree/master/Code/MinimalLib#live-demos)


## Future plans