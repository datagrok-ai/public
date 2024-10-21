
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
  scalars, help, history. Dock manager is used as a container,
  providing layout customization.

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

Additional RFV view variations with Model Hub integration will be
available (TODO):

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
  `Reactive Tree Driver`, `Rich Function View` and `Tree Wizard`.

- `Compute` package is the primary consumer of `compute-utils`
  library. it depends on `compute-utils`.

- `compute-api` is a library which provides a way for other packages
  to depend on a dynamically loaded version of `compute-utils`. The
  code is loaded in the `compute` package initialization, so for it to
  work, `compute` must be loaded before using `compute-api`. It is
  used to update `Standalone Function View` without a need to rebuild
  and redeploy custom JS UI models.

- `webcomponents` library has wrappers around platform forms and
  viewers. It also includes a web component version of docker
  WM. (TODO: move docker WM to a separate package, or even merge with
  the platform core one).

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

2. `WebComponents` initialization is done during `Compute` init.

3. `compute-api` initialization is done during `Compute` init.

4. `Compute` initialization is done. Model Hub integrated abstractions
   can be launched, since all dynamic dependencies are ready.


Note that `Compute` depends only on `webomponents-vue` library, which
doesn't include web components themselves, that is why `WebComponents`
init is necessary. Also a global platform vue.js library must be used.
