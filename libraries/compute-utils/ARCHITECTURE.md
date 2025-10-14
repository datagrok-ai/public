
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
  also includes fitting and optimization related code with ui.

- `Compute` Provides `Model Hub` app via package functions.

- `Compute2` contains Vue.js UI code for, `Rich Function View` and
  `Tree Wizard` and provides Vue.js based editors.

- `compute-api` is a library which provides a way for other packages
  to depend on a dynamically loaded version of `compute-utils`. The
  code is loaded in the `Compute` package initialization, so for it to
  work, `Compute` must be loaded before using `compute-api`. It is
  used to update `Standalone Function View` without a need to rebuild
  and redeploy custom JS UI models.

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

1. `Compute` package is starting. This could be triggered via an URL
   or via starting Model Hub app from UI.

2. `compute-api` initialization is done during `Compute` init.

3. `Compute` initialization is done. `compute-api` dependent code
   can be launched, since all dynamic dependencies are ready.

4. `Compute2` initialization is done when starting a `Rich Function
   View` or `Tree Wizard`.

5. `WebComponents` initialization is done during `Compute2` init.

6. `Compute2` initialization is done, `Rich Function View` or `Tree
   Wizard` based views can be rendered.


Note that `Compute2` depends only on `webomponents-vue` library, which
doesn't include web components themselves, that is why `WebComponents`
init is necessary. Also a global platform shared Vue.js library is
used.
