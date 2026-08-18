# U2 changelog

## v.next

* Tokens: added semantic selection (--dg-selection-bg/text), elevation scale (--dg-shadow-1/2/3), focus ring (--dg-focus-ring), and overlay radius (--dg-radius-l) tokens; component CSS converted
* Tokens: adopted the P7 designer-pass values as defaults — layered elevation shadows, 4px/6px radii, accent-tinted selection, 2px focus ring, calibrated success/failure colors, 16/14px heading hierarchy
* Added theme-dark.css (parked draft, [data-dg-theme='dark']); gallery gained a theme toggle

* Added item actions (list-item-rendering recipe): the Action type, rowActions() hover/focus-revealed icon block (opacity-hidden so buttons keep their space and stay tabbable), actionsMenu(), and VirtualList contextActions — right-click selects the row and opens the full action list at the cursor
* Added timestamp() element factory: locale short date (year only when not current), full date-time in the tooltip; accepts Date, numbers, strings, or anything with toDate()
* Added css/adaptive.css: u2-adaptive container + u2-p1/u2-p2 priority classes — secondary row content hides as the pane narrows (container queries, fixed 420/340px defaults)
* Menu: item icons now render as icon glyphs instead of raw text

* Overlay: fixed positioning strategy + animationFrame autoUpdate (dock_spawn moves anchors outside scroll/resize); Combobox popup anchored to the input, not the stretched root
* host(): optional closeIn dock manager — auto-detaches on pane ✕ via onClosed (the platform never kills non-view dock content)
* Added ambient scope ownership (Scope.ambient/runWith, Component.build) — components and signal bindings created under an owner are disposed with it
* Added fluent element factories (div/divV/divH/span/label/h1-h3/link/button/bigButton/panel/empty) accepting T | Signal<T>
* Added Input<T> foundation (label/editor layout, validators + validity signal, inline variant) with TextInput/TextArea/BoolInput (+switch)
* Added Tooltip service: singleton in the shared overlay layer, 300ms delay with instant warm switch
* Added Splitter: VS Code-style sashes, pointer-capture drag, keyboard resize, nesting
* Added TabStrip: lazy panes, dirty markers, middle-click close, overflow dropdown
* Added Menu: fluent builder (DG.Menu-shaped), submenus with hover intent, typeahead, machine-state signal
* Added VirtualTree over VirtualList: lazy branches, inline rename, expandPath, 100k-node scale
* Gallery: per-page scope — navigating away disposes the previous page's components
* Added ChoiceInput (styled native select, platform convention) and MultiChoiceInput
* Added NumberInput: int/float parsing, min/max clamp on commit, hover spinner, arrow/wheel stepping
* Added Dialog: fluent DG.Dialog-shaped API, modal backdrop, focus trap, dragging, stacking, per-show scope
* Added Accordion (lazy panes on per-pane scopes) and Breadcrumbs (overflow ellipsis collapse)
* Added Toolbar: flat buttons/toggles/menus/separators/slots, overflow-collapse into a mirror Menu
* Added Form: aligned label column, aggregate validity signal, validate/getValues/setValues, condensed/wide
* Overlay: shared monotonic layer counter — later-opened popups/dialogs always paint on top (fixes combobox-in-dialog occlusion)
* Added PropertyGrid: two-column inline editors from PropDescriptors, collapsible categories, write-through values signal
* Added AsyncView + loader/skeleton: standardized loading/empty/error-with-retry rendering over AsyncSource
* Added Perf instrumentation and the gallery performance page (budgets: 200-widget construct <15ms, 1M-row 60fps scroll, zero-leak round-trip — all passing)
* Added node test harness: tests/dom-shim.js + 25-test smoke suite (npm test), leak assertion on every test
* Added dg input bridge (asDartInput → JsInputBase) and columnInput/tableInput pickers
* VirtualList: rowRole + render row argument + keyOf identity-preserving selection; VirtualTree native per-row ARIA (retrofit removed)
* Menu: Home/End navigation, focus restore on close; Input: name option (Form keys values by it)
* Added dg-ui/1 spec layer: Registry + manifest, renderSpec/SpecContext with $.path bindings and cmd: commands (no eval), SpecInstance.dump round-trip, node-level error containment; 13 components registered; manifest.json + npm run manifest; gallery spec page
* Spec: adopt()/createWithChildren child hooks + childProps metadata + 'json' prop type — splitter/accordion/tabs/form now spec-registrable (manifest: 16 components)
* Added schema-driven forms (dg): propertyForm/objectForm over PropertyLike — editors by type, required from nullable, min/max, overrides; ObjectForm.refresh()
* Overlay.nextLayer(): u2-local monotonic order inside the overlay host's own z band (platform surfaces keep their banded z-orders)
* host() stamps data-kill-on-close — the new core dock contract kills the component on pane ✕ (older cores covered by the closeIn onClosed auto-wire)
* Added fromDartInput (dg): any DG.InputBase participates in u2 Forms — value sync, Dart validators feeding Form.validity via the new core validity surface (feature-detected), full platform chrome; asDartInput/fromDartInput unwrap instead of double-wrapping
* objectForm: editors 'auto' mode — real DG.Property fields resolve the platform's own editor (DG.InputBase.forProperty) wrapped via fromDartInput, generated u2 editors as fallback
* Added Wizard: step rail with visited/current states, canProceed gates composing Form.validity (reason shown), lazy per-step scopes, Enter-to-advance, openInDialog
* Added TypeAhead<T>: complex rendered option rows (per-render child scopes), itemText filtering, selected signal cleared on typing, async loading/empty/error states, minChars; dg userInput() over dapi.users with avatar rows
* Added badges: badge (5 token variants), countBadge (99+ cap, hides at 0, signal-driven), status dot, removable tag
* Added ProgressBar: clamped 0..1 value signal, runtime-switchable indeterminate, percent label, reduced-motion fallback
* Added BasicTable<T> + tableFromMap: semantic compact table for small data, string/element cell renderers under per-render scopes, live signal rows, keyboard row selection
* Added MultiSelect: popup multi-select Input<string[]> — tag field with +N overflow, stay-open checkbox popup, auto filter, select all/none over visible rows
* Added icons: icon(name) emitting the platform's grok-icon fal fa-* convention (FA Pro in-app); vendored FA Free 5.15.4 woff2 + .u2-standalone-gated codepoint shim for standalone hosts (1002 solid names)
* Added iconButton (flat square, toggle signal) and buttonWithIcon (optically aligned to the button typography)
* Added ButtonGroup: segmented/flat, toggle none|single(radio)|multi, iconOnly with auto-tooltips, form|toolbar|ribbon densities, roving tab stop tracked by id
* Icons: standalone fallback upgraded to the platform's own FA Pro 5.15.4 (true Light face, full 2309-name map, brands included); faces requested on first icon() use to avoid blank first paint
* Added dropDownButton: plain (whole button opens a fresh Menu per open) and split (main action + attached arrow) modes, aria-expanded via a self-disposing per-open scope, ArrowDown opens
* Added MenuBar: classic top-menu over Menu popups — click opens, instant hover-switch while open, roving arrows switch the open menu, transparent chrome for ribbon/toolbar hosting
* Added DateInput + DateTimeInput (the 18th audited control): shared APG-pattern calendar popup, strict yyyy-MM-dd[ HH:mm] parse/commit with min/max clamp, native Date values, time segments keep time-of-day on day change; objectForm datetime properties now generate DateTimeInput (Date/dayjs/ISO coerced)
* Added AsyncPager<T> (src/core/pager.ts): the paged twin of AsyncSource — pages accumulate into one items signal, short page marks done, reset() drops in-flight responses and re-reads the total; dg dapiPager() pages any grok.dapi.* collection (factory-per-request, thunk filter, order pass-through); gallery pager page

* Initial scaffold: library skeleton, platform-free-core eslint guard, vendor/ layout, standalone gallery shell
* Added design tokens (css/tokens.css): 74 --dg-* tokens codifying current Datagrok styles, compact density axis, gallery token sheet
* Added signals core: vendored @preact/signals-core 1.14.4, Scope/Component lifecycle with leak counter, bindText/bindValue with echo suppression, dg adapter (host(), toSignal/toObservable, leakReport)
* Added VirtualList (src/components/list.ts): VS Code listView virtualization ported onto the u2 core — fixed row height, row recycling, signal-driven selection, keyboard navigation, listbox ARIA; css/list.css and a 1M-row gallery page
* Added AsyncSource (src/core/async-source.ts): debounced, abortable idle/loading/ready/empty/error-with-retry contract for every async control
* Added Overlay (src/core/overlay.ts): single body-level layer host, @floating-ui/dom positioning (flip + size), dismissal on outside pointerdown / Esc / anchor detach
* Added Combobox (src/components/combobox.ts): hand-written WAI-ARIA 1.2 combobox machine (Zag.js vendoring evaluated and rejected — see the file header), sync and async sources, css/combobox.css and a gallery page
