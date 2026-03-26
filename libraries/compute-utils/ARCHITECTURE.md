
This document describes a high level design of compute related
components and theirs primary goals.


## Types of FuncCall based abstractions

The primary goal of `compute-utils` and related packages/libraries is
to provide abstractions on top of FuncCalls.

- `Model Hub` - a list of scripts with a special annotation which are
  grouped by several nested categories. It is the entry point for end
  users.

- `Rich Function View` - is a Vue.js framework component which
  generates UI for FuncCalls. It includes input form, output viewers,
  scalars, help, history. Dock manager is used as a container for
  components. Also provides integrations with fitting and optimization
  tools.

- `Reactive Tree Driver` - is a configurable state store, backed by
  FuncCalls, and a linking mechanism between those states. For an end
  user it will be a sequence of steps workflow, with some configurable
  runtime mutations capabilities. It is created via a typescript
  configuration with simple link handler functions and provides many
  addition features: runtime tree structure mutations, tree queries
  for handlers, action handlers, different types of links (validators,
  metadata, tree mutations), consistency checking, readonly subtrees,
  full and partial saving/loading. For the UI integration driver
  provides (using rxjs observables) main state tree, separate meta
  states for tree nodes, and a message api for tree mutations.

- `Tree Wizard` - is an UI layer for `Reactive Tree Driver`, each
  FuncCall is rendered the same as in `Rich Function View`. `Tree
  Wizard` has a navigation tree and additional driver integration
  logic/ui elements. Also integrates into Model Hub.

- `Standalone Rich Function View` - an integration of `Rich Function
  View` with Model Hub. Provides a simple way to make FuncCall UI only
  using annotations.

- `Standalone Function View` - an integration of history only `Rich
  Function View` with Model Hub. Provides a way to make completely
  custom FuncCall UI, still integrated with Model Hub.



## Libraries and packages overview

The implementation of mentioned above abstractions is split across
multiple packages and libraries.

- `compute-utils` library is the main source of all code related to
  dealing with FuncCalls. The main things here are `Model Hub`,
  `Reactive Tree Driver`, `Rich Function View` and `Tree Wizard`. It
  includes fitting and optimization related code with ui. Also
  includes leagacy views, provided by `Compute` package.

- `Compute` provides legacy function views for old packages support.

- `Compute2` contains Vue.js UI code for, `Rich Function View` and
  `Tree Wizard` and provides Vue.js based editors. It also provides
  `Model Hub` app via package functions.

- `compute-api` is a library which provides a way for other packages
  to depend on a dynamically loaded code of `compute-utils` in
  `Compute` pacakge, it is mostly used for legacy `Compute` views. For
  `Compute2` it provides typings and `Standalone Function View`.

- `webcomponents` library has wrappers around platform forms and
  viewers.

- `dock-spawn-dg` provides dock manager web component.

- `webomponents-vue` library is a consumer of `webcomponents` library
  typings, providing typed vue.js components.

- `WebComponents` is a package which consumes `webcomponents` and
  provides webcomponets in the DOM. A special initialization function
  must be called by a consumer package before using any of web
  components.

- `LibTests` contains an extensive `Reactive Tree Driver` logic test
  suite (100+ tests).

Here is an initialization sequence:

1. `Compute2` package is starting. This could be triggered via an URL
   or via starting Model Hub app from UI. It will also make avaliable
   `Standalone Function View`.

2. `Compute2` will initilize `WebComponents` package, making
   webcomponents avaliable in DOM. If present `Compute` will be
   initilized as well.

3. `compute-api` legacy views code initialization is done during
   `Compute` init, and `Standalone Function View` code initialization
   is done in `Compute2` init.

4. The model code and package are loaded and initilized, model related
   FuncCalls are created.

5. Initialization are done, `Rich Function View` or `Tree Wizard` can
   be rendered.


Note that `Compute2` depends only on `webomponents-vue` library, which
doesn't include web components themselves, that is why `WebComponents`
init is necessary. Also a global platform shared Vue.js library is
used.
