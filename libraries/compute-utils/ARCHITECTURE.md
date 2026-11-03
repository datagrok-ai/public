# Compute components architecture

This document describes a high level design of compute related
components and their primary goals.


## Types of FuncCall based abstractions

The primary goal of `compute-utils` and related packages/libraries is
to provide abstractions on top of FuncCalls.

- `Model Hub` - a list of scripts with a special annotation which are
  grouped by several nested categories. It is the entry point for end
  users.

- `Rich Function View` - is a Vue.js framework component which
  generates UI for FuncCalls. It includes input form, output viewers,
  scalars, help, saved runs (Runs panel) and run comparison. Dock
  manager is used as a container for components. Also provides
  integrations with fitting and optimization tools.

- `Reactive Tree Driver` - is a configurable state store, backed by
  FuncCalls, and a linking mechanism between those states. For an end
  user it will be a sequence of steps workflow, with some configurable
  runtime mutations capabilities. It is created via a typescript
  configuration with simple link handler functions and provides many
  additional features: runtime tree structure mutations, tree queries
  for handlers, action handlers, different types of links (data,
  validators, pipeline validators, check links, rules, view metadata,
  node metadata/selectors, tree mutations, funccall actions, return
  links), consistency checking, readonly subtrees, full and partial
  saving/loading, headless run and save without UI. For the UI
  integration driver provides (using rxjs observables) main state
  tree, separate meta states for tree nodes, and a message api for
  tree mutations.

- `Tree Wizard` - is a UI layer for `Reactive Tree Driver`, each
  FuncCall is rendered the same as in `Rich Function View`. `Tree
  Wizard` has a navigation tree and additional driver integration
  logic/ui elements. Also integrates into Model Hub.

- `Standalone Rich Function View` - an integration of `Rich Function
  View` with Model Hub. Provides a simple way to make FuncCall UI only
  using annotations.

- `Standalone Function View` - a completely custom FuncCall UI
  (`CustomFunctionView`), still integrated with Model Hub. Its Save and
  Open dialogs come from the host's run storage (Compute2 provides
  them with My runs and spaces).



## Libraries and packages overview

The implementation of mentioned above abstractions is split across
multiple packages and libraries.

- `compute-utils` library is the main source of all code related to
  dealing with FuncCalls. It contains `Reactive Tree Driver`,
  `Model Hub` (`model-catalog/`), run storage and spaces
  (`history-utils/`), and `function-views/` with fitting, sensitivity
  analysis and optimization views and `CustomFunctionView`. Legacy
  views and components live in `old-views/` and `old-components/`,
  used only by `Compute` and `compute-api` consumers.

- `Compute` provides legacy function views for old packages support.

- `Compute2` contains Vue.js UI code for `Rich Function View` and
  `Tree Wizard` and provides Vue.js based editors. It also provides
  `Model Hub` app, run storage UI (Save run dialog, Browse Runs), run
  comparison (Compare Runs) and headless `RunWorkflow` /
  `SaveWorkflowRun` functions via package functions.

- `compute-api` is the public contract for other packages. It
  provides factories for legacy views from dynamically loaded
  `Compute` code (`window.compute`), `Reactive Tree Driver` types and
  state guards, a headless workflow API (`startWorkflow`,
  `runWorkflow`, `saveWorkflowRun`, calling `Compute2` functions) and
  the `CustomFunctionView` base class.

- `webcomponents` library has wrappers around platform forms and
  viewers.

- `dock-spawn-dg` provides dock manager web component. It is
  registered as `<dock-spawn-ts>` by the `WebComponents` package, so
  `WebComponents` init is needed for it as well.

- `webcomponents-vue` library is a consumer of `webcomponents` library
  typings, providing typed vue.js components.

- `WebComponents` is a package which consumes `webcomponents` and
  provides web components in the DOM. A special initialization function
  must be called by a consumer package before using any of web
  components.

- `LibTests` contains an extensive `Reactive Tree Driver` logic test
  suite.

Here is an initialization sequence:

1. `Compute2` package is starting. This could be triggered via a URL,
   via starting Model Hub app from UI, or by calling any `Compute2`
   function: editors, Browse/Compare Runs menu items, or
   `startWorkflow` / `runWorkflow` from `compute-api`.

2. `Compute2` will initialize `WebComponents` package, making
   web components available in DOM. If present `Compute` will be
   initialized as well.

3. Calling `Compute:init` (directly or via `initComputeApi()`) loads
   the `Compute` bundle, which exposes legacy views code as
   `window.compute` for `compute-api`.

4. The model code and package are loaded and initialized, model related
   FuncCalls are created.

5. Initialization is done, `Rich Function View` or `Tree Wizard` can
   be rendered.


Note that for web components `Compute2` depends only on the
`webcomponents-vue` library, which provides typed wrappers but doesn't
include or register the web components themselves, that is why
`WebComponents` init is necessary. Also a global platform shared
Vue.js library is used.
