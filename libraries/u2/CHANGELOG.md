# U2 changelog

## v.next

* Testing: every named spec node stamps `data-u2-name` on its root element — one automation id shared by browser checks, package tests, tutorials and agents — and `e2e/` runs the designer checklist against the serverless client (`?mode=local`, no login, no publish) in ~25 s via `npm run e2e:local`, with the stand-only checks (publish, URL route, Browse card) kept apart in `npm run e2e:stand`
* Designer: the property panel heads itself with the node's name, tag and — for a node that failed to build — that it is broken; a prop the component keeps neither as a signal, a member nor an `options` entry (a `create` wrapping a bare button, a splitter storing `direction` as `_horizontal`) now reports the value the spec set instead of a blank; and selecting the row that is already selected re-asserts the shell object, so a context panel that missed an update comes back on a second click
* Designer: a read-only inspector over any rendered spec — `designerView(spec)` is a platform view whose canvas is the live `SpecInstance`, whose toolbox is its structure tree, whose ribbon carries Open…/Copy spec and the Design/Run toggle, and whose status bar shows the selection path, the node count and how many nodes are broken
* Designer: selection is the platform's own mechanism — clicking a node (on the canvas or in the tree) makes a `SpecNodeRef` the shell's current object, and `registerSpecNodeHandler()`'s `ObjectHandler` renders it in the context panel: identity, the component's properties grouped by category with their live values, every event the component declares with its `cmd:` wiring, and the node's bindings verbatim; a node that failed to build is selectable too, and the panel says why
* Designer: Design mode is one transparent glass pane over the canvas — components render exactly as they run but take no clicks, and the pane hit-tests through `SpecInstance.nodeAt` to outline the node under the pointer and the selected one; Run mode removes the pane and the outlines, same rendered tree either way
* Property metadata: `PropertyLike` gained `category`, which the property panel groups by — it reaches the panel through metadata-generated `getProperties()` with no per-component code
* dg: `appView({toolbox})` puts a component or an element in the view's toolbox pane, alongside the existing ribbon and status-bar options
* **Breaking** — Component/Control split: `Component` is now the DOM-free base (name, scope, effects, disposal) and `Control extends Component` is everything visual, so a data source can be a component on a form without pretending to render; every u2 component is a `Control`, and `Component.build` is `Control.build`
* Spec: a node can carry a `name` (unique within the spec, preserved by `dump()`, a duplicate warns instead of failing), and `SpecInstance` exposes what it rendered — `nodes()` over every node including plain HTML ones, `node(name)`, and `nodeAt(element)`, the hit-test a designer selects through
* Widget convergence: a component answers the questions the platform asks every widget — `getProperties()` generated from its registry metadata over its own signals, `getFunctions()`, `onEvent()` (component events plus the DOM events of that name), `aiDescription` and `getWidgetStatus()` — through the platform-free `WidgetLike`; `host()` now returns a real `DG.Widget` subclass delegating all five to the component, with its functions minted as platform `Func`s, and a compile-time assertion keeps the structural copy from drifting off `DG.Widget`
* Spec registry: a component prop is described the way the platform describes every property — `PropertyLike` (the `DG.Property` shape, now in the platform-free core) plus `bindable`/`twoWay` — and prop types are platform TYPE strings: `float` → `double`, `string[]` → `string_list`, `json` → `object`
* **Breaking** — PropertyGrid: a row descriptor's `type` speaks that same platform vocabulary (`float` → `double`), so the manifest describes one type vocabulary throughout
* Widget convergence: `getProperties()` reports what a control actually holds — a same-named signal, else a same-named member, else the matching `options` entry — instead of `undefined` for every prop whose state lives in options; writes stay signal-only, so a prop without one renders read-only
* Fixed: a component-event subscription left open survived `dispose()` and kept delivering (the DOM half already died with the scope)
* Fixed: Designer — Open… left an arbitrary node of the new spec selected (and announced it) whenever a row had been selected before: the tree preserves selection by row key, and ordinal ids collide across specs; the tree gained `VirtualTree.clearSelection()`, and `open()` clears before it swaps roots, so the new root is what ends up — and stays — selected

* Inputs: a validation message is a line of its own under the editor, always — beside a 100px color field it used to fit on the row and render off the edge of a narrow pane, so the same message appeared under one input and out of sight next to another
* Inputs: commitOn — an input can report its value when the edit is committed (blur or Enter), as Dart's fireChanged does, instead of on every keystroke; asked for where the consumer acts on the value irreversibly
* dg: propertyEditor() renames on commit — the name field reports when it is left or Enter is pressed, so a caller keying properties by name no longer renames the property through every prefix of what is being typed (`dose` → `n` → `na` → `nam` → `name`)
* Tooltip: a hovered element torn down before the hover delay elapses takes its pending tooltip with it, and a detached anchor is never shown against — a rebuilt list used to leave the tooltip of a row that no longer exists stuck over the page
* Inputs: every input now carries the platform options rail — units and side icons sit right of the editor on one 28px row whose underline continues the field's, with the units pinned leftmost (Dart InputBase.options); the rail travels with the editor through the platform bridge, so a u2 value editor keeps its units and icons inside a Dart form
* NumberInput: the slider and the − + clicker are the platform chrome exactly — hidden at rest, revealed on hover anywhere over the input, the clicker 20px and blue on the rail before the units, the slider spanning the field width under it (overlapping the next row by 9px, as d4.css positions it); u2's own hover spinner has no platform counterpart and is now opt-in (spinner: true), so a bare number input is a bare field
* SliderInput: dropped the readout beside the track — a bare track is what the platform shows; the value surfaces in the tooltip, on hover and while dragging
* TextInput: the search variant's clear icon shows only while the pointer is over the input, at the right edge of the field, as the platform's does
* DateInput/DateTimeInput: the calendar affordance is the platform's calendar-alt icon on the options rail, always visible, beside a 140px field
* ColorInput: the platform layout — a 100px hex field with the swatch on the options rail right of it; swatchOnly drops the field and the rail underline with it
* dg: the columnsInput check-list is an anchored popup instead of a modal dialog — it opens under the control (up, near the viewport bottom) as Dart's own selector does, and follows its dismissal: an outside click or Esc cancels, OK or Enter commits (from the search box too, so arrow → Space → Enter picks a column without touching the mouse), and clicking the control closes it again; the list is a whole number of rows tall, so it never ends on a half-drawn one
* dg: propertyForm/objectForm route on the property's own editor hints, as the platform does — inputType first (Radio, Color, Font, Image, Tags, Slider, Choice, MultiChoice, TextArea, Switch…), then editor (textarea, password, switch, slider), then the type; list, map and file properties get the list/multi-choice, map and file editors instead of a disabled text box
* dg: a u2 value editor no longer repeats the units of the property the platform bound it to — the Dart host renders those itself, so the bridge drops the rail postfix (and the rail, where nothing else is on it) on that path alone
* dg: the options rail is pinned to the platform's own size instead of `font-size: smaller`, which resolved differently under a Dart form's root than in a standalone page
* dg: rsaInput() accepts a dragged key whose `connection` cannot be read — js-api's `FileInfo.connection` throws on deployed bundles, which silently rejected every drop — and says so when the file itself cannot be read
* dg: propertyEditor() captions every field in sentence case and takes per-field validators, so a caller can report what only it knows (a name already taken) on the field instead of dropping the edit
* dg: propertyEditor() — a reusable editor of property METADATA (the `IProperty` options themselves), with the platform's field vocabulary and its applicableTo visibility rules: an identity section plus a type-dependent one that swaps min/max/step/sliders for choices as the type changes; every edit reports the whole record, which the caller turns back into a property with DG.Property.fromOptions
* dg: tableInput() and tablesInput() gained the platform input's import action — a folder-open icon on the options rail opens a local file, adds the table to the workspace and picks (or checks) it; csv/tsv/txt parse and d42 deserializes, and an extension neither reader takes is refused with the platform's own message, since js-api has no door to the FileHandler that opens xlsx, sdf and the rest for the Dart input
* dg: rsaInput() takes the dragged object however the platform hands it over (a DragDropArgs or the object itself) and leaves a drag that is not a key file of a file share alone, instead of reading whatever was dropped
* dg: the file and files pickers sit on the shared options rail like every other input's side icons, so their underline continues the field's and they cross the platform bridge with the editor; the file, files and rsa rails carry the platform's icon metrics — a blue 16px icon box, 8px from the one before it, in a control-height row
* TypeAhead: a typed query highlights its first match, so plain Enter accepts the suggestion on screen instead of doing nothing; TagsInput still highlights nothing (its box's own text is what Enter is for)
* Inputs: a validation message with no room beside the editor takes the next line whole instead of shrinking to one word per line
* dg: the name-valued column pickers — columnInput (both variants), columnsInput, columnsMapInput and aggregatedColumnsInput — now follow onColumnNameChanged and REMAP the name they hold; a rename fires that event alone, so until now every one of them silently pointed at a column that no longer answered to it
* MultiChoiceInput: emptyText — an empty list reads as empty instead of as an empty box; tablesInput says 'No open tables'
* dg: the columns dialog gained All / None links and an 'N selected' count, and its checkboxes name their column for a screen reader
* propertyForm/objectForm: a datetime field no longer reports a change (and no longer writes the property back) the moment it is constructed — every read builds a fresh Date, which identity comparison took for an edit
* dg: fileInput() says 'File does not exist.' without repeating the path the box already shows
* NumberInput: the inline slider draws a focus ring when it is reached by keyboard, as the standalone slider does
* Inputs: a computed-style sweep over the side-by-side demo closed the measured gaps against the platform's own metrics — the choice box at 200px with its 20px dropdown inset, the text box at 100px, bigint and qnum at 140px, radio and multi-choice stacked as 20px rows with the box 4px in and 8px off its label, map rows on the platform's key/value/button table, the font row at its own 11px scale, and one 16px blue icon box for every glyph on the options rail
* Inputs: the platform draws its side affordances with icons, so the number clicker, the map row buttons and the font size steppers do too, instead of text plus and minus signs
* FontInput: the − + size steppers the platform puts beside the size box
* ListInput: the expand toggle moved onto the options rail, where the platform's is, and the field stopped drawing a second underline under the editor's
* ImageInput: the preview is the platform's fixed 300×200 box, with REMOVE as a link over its top-right corner
* TextInput: the password reveal icon sits over the right edge of the field, as the platform's does, instead of taking width from it
* dg: an editor hint applies only to the types the platform honours it for — `{editor: 'textarea'}` on an int property builds the int editor, as `InputBase.forProperty` does

* dg: the input bridge fires the platform events from a microtask instead of from the value effect — a platform onChanged subscriber that writes a u2 signal back no longer trips the signal runtime's cycle guard and kills the bridge; one edit is still one fireInput/fireChanged pair
* dg: the bridge re-registers the u2 validity feed on the dart handle the platform re-homes the input on, so a registered value editor's invalid state now marks the platform input invalid and disables a dialog's OK
* Inputs: Input.system() marks a value the input prunes on its own (an item list that lost the selected entry, a column dropped from the table) as programmatic, so the bridge reports onChanged alone, as Dart does
* dg: inputForProperty() takes assumeWritable — a value editor gets the editor the type asks for even though a FuncParam carries no setter (a form still renders a setter-less property read-only)
* dg: fileInput() — a value written from outside drops the verdict on the path it replaces and cannot be overruled by a check still in flight; a path the listing does not carry now reports "does not exist" instead of staying silently valid; the extension filter is applied in file mode only, so a folder whose name carries a dot resolves
* dg: filesInput() and fileInput() cancel a read still running when the input is disposed, instead of settling onto a disposed component
* dg: rsaInput() refuses a drop of several files, as the single-file input does, and ignores a platform drop that lands after dispose
* dg: aggregatedColumnsInput() rows are built from ChoiceInput and icon buttons instead of hand-rolled select/button markup

* dg: added tablesInput() — a checkbox per open table with the names as the value, following table open and close (port of Dart TablesMultiChoiceInput)
* dg: columnInput() gained the rich variant — a type-ahead over the columns with type glyphs, names and semantic types (Dart's ColumnComboBox affordances), on the same name-valued contract as the plain choice variant
* dg: added columnRenderer() — an ObjectRenderer over column names any u2 control can take; a name the table no longer has still renders
* dg: added columnsInput() — Dart's `(3) age, sex, race` / `(N) All` summary opening a searchable check-list (u2 VirtualList, not the platform selectColumns modal) with the picked columns on top; OK commits, cancel does not, changeTable() starts over (port of Dart ColumnsInput minus its additional-column checks)
* dg: added columnsMapInput() — a column picker per key, the key's type restricting what it offers (int and double count as one type, as in Dart); unmapped keys stay out of the value (port of Dart ColumnsMapInput)
* dg: added aggregatedColumnsInput() and aggregationsFor(type) — column + aggregation rows whose aggregation list is ddt's per-type registry minus the exclusions of Dart's aggregation editor, so string and bool columns offer the counts and datetime those plus range (port of Dart editColumnAggregations)
* MultiChoiceInput: setItems() — the checked items that are still there stay checked, as in the Dart input
* dg: added fileInput() — path box with a local-file picker and a drop zone, checked against the file share as it is typed (debounced, one check at a time, late answers to abandoned paths dropped) with the platform's empty / still-checking / does-not-exist states; a connection namespaces the path with its nqName (port of Dart FileInput)
* dg: added filesInput() — several files at once, each row with its own percentage and a ✕ that cancels the read or drops the file, re-adding a name replaces it, and the loaded FileInfos are published once every read has settled (port of Dart FilesInput)
* dg: added rsaInput() — a key file goes in and is held masked: its text when it opens with a PEM header, BASE64 of its bytes when it does not; keys dragged from the platform's file browser are read through grok.dapi.files and must be PEM (port of Dart RsaInput)
* Added FontInput: the platform's strict font string ("<weight> <style> <size>px <family>") through a size box with a preset popup, bold/italic toggles and a family select; an unknown family survives as a hidden option, generic families stay lowercase and unquoted (port of Dart FontInput)
* Added ImageInput: preview of the URL or data URL it holds, with a REMOVE link revealed on hover (port of Dart ImageInput)
* TextInput: password variant (masked field plus an eye toggle) and autoResize (width follows the text between minWidth and maxWidth); the platform PASSWORD_PLACEHOLDER focus/blur dance stays out of the core input
* Added TagsInput<T>: object tags over AsyncSource and ObjectRenderer — chips with per-item deleteEnabled, sync or async suggestions that drop what is already picked, allowNew creation, and the full keyboard contract (arrows/Home/End, Enter accepts or creates, Esc closes the popup only, Backspace pops the last chip); MultiSelect stays the string-list control
* TagsInput: Ctrl+Enter (and Cmd+Enter) no longer adds a chip and no longer swallows the key — it stays the hosting dialog's OK shortcut, as in the Dart input
* FontInput: the size popup toggles — a click on the box it is already under closes it, and the next click brings it back
* Overlay: closing an overlay now drops its cleanup from the owning scope, so an input reopened all day stops accumulating dead closures (Scope.disown)
* Added SuggestionList (src/core/suggestion-list.ts): the overlay-hosted listbox, the active-row machine and the loading/empty/error rows now live once and back both TypeAhead and TagsInput
* Added QNum (src/core/qnum.ts): the qualified-number packing ported from ddt's QNum — create/getValue/getQ/parse/format over the two lowest mantissa bits, byte-for-byte with the Dart original
* Added QNumInput: editor over a packed qnum double, taking a leading < or >; unparseable text stays and marks the input invalid (Dart reports null for it)
* Added BigIntInput: integer of any size kept as text end to end, so no digit ever passes through a double (port of Dart BigIntInput)
* Added DynamicInput: a read-only input that renders its value through a supplied renderer and swaps the element on every write (port of Dart DynamicInput)
* dg: added metaInput() — DynamicInput over the object handler's card, with Dart MetaInput's `<null>` marker
* NumberInput: platform parity — min/max now validate instead of clamping typed text (the typed number reaches the value signal and the message shows; only the slider and the stepper stay inside the bounds), plus the showSlider/showPlusMinus equivalents (slider/clicker), a display format hook, and Dart's dynamic float precision
* Inputs: postfix option on every input — units after the editor, as InputBase.addPostfix
* propertyForm/objectForm: number properties now carry step, format, units, showSlider and showPlusMinus into the generated editor with the platform's defaulting rules; inputForProperty() exposes that mapping
* propertyForm/objectForm: bigint and qnum properties now generate BigIntInput and QNumInput instead of a number input
* dg: added dartInputFor(build) — a value editor that builds its u2 input once the platform has bound the property, so registered editors keep the property's metadata
* Added SliderInput: bare range track, (max - min) / 100 default step (port of Dart SliderInput)
* Added RadioInput: native radio group sharing one name, setItems, per-item tooltips, raised-button variant (port of Dart RadioInput)
* Added ColorInput: hex field plus a swatch fronting the native color picker, swatchOnly variant; bad hex marks invalid and keeps the last good color (port of Dart ColorInput, platform palette not ported)
* Added ListInput: comma-separated list with Dart's quote-aware split and the expand-to-textarea toggle
* Added MapInput: key/value rows with per-row add/remove; duplicate keys go through the validity signal and empty-key rows stay out of the value
* MultiSelect: allowCustom — Enter in the filter box adds the typed text as a tag, matching existing items case-insensitively
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
