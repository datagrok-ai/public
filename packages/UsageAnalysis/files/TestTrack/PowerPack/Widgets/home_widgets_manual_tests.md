---
feature: widgets
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [powerpack.cp.welcome-view-renders, powerpack.cp.spotlight-widget-renders, welcome-view-power-search]
realizes: []
realized_as:
  - home-widgets.test.ts
related_bugs: []
---

# Home page Widgets — manual test cases

Coverage of the **widgets on the Home page** (Datagrok landing view): rendering, management (close/Customize),
persistence, search, navigation, and the functionality of individual widgets (Spotlight, Community, Usage/Reports).

Context, environment, selectors, and pitfalls are in `widgets_tests_context.md`.
Cases are consolidated and structured for a direct port to Playwright (one case = one end-to-end flow).

---

## 1. Environment and common preconditions

- **Instance:** `https://dev.datagrok.ai/` under a logged-in user (`opavlenko+playwright`, has the admin role).
- **Home Page** is open (via login or the **Home** icon). The **Spotlight**, **Community**, **Usage**, **Reports** widgets are visible.
- The main end-to-end check for all cases is **no errors** (`pageerror`, `console.error`, error balloons).
- **Cleanup:** cases that change widget visibility restore the initial state (all widgets `ignored:false`).

---

## 2. Display — rendering the Home page and widgets

### Widgets-Display-01 — Home page renders with all widgets
**Type:** smoke
**Preconditions:** section 1.
**Steps:**
1. Log in / open the Home page → reaction: the `.power-pack-welcome-view` view opens, `.power-pack-widgets-host` appears with child `.power-pack-widget-host` elements.
2. Inspect the rendered widgets (presence, order, titles, content).
**Final check:**
- ✅ a search box is present (placeholder starts with "Search everywhere");
- ✅ the **Spotlight**, **Community**, **Usage**, **Reports** widgets are visible (`[widget-title="..."]`);
- ✅ **order:** Spotlight (order -1) comes before Community (order 6) in DOM order;
- ✅ **titles:** Community shows `.d4-dialog-title` = "Community"; Spotlight has no title (`showName=false`, `.d4-dialog-header-hidden`);
- ✅ each widget's content (`.power-pack-widget-content`) is non-empty;
- ✅ no errors (a crashed widget logs `console.error` and removes its host — this must not happen).
**Postconditions:** —

---

## 3. Controls — managing a widget (hover + close)

### Widgets-Controls-01 — Hover shows the close icon, click removes the widget
**Type:** functional
**Preconditions:** Home page open, Community visible.
**Steps:**
1. Hover over the **Community** widget → reaction: the close icon `i.grok-font-icon-close` becomes visible (tooltip "Remove").
2. Click the close icon → reaction: the widget disappears from `.power-pack-widgets-host`.
**Final check:**
- ✅ on hover the close icon is visible, without hover it is hidden;
- ✅ after the click `[widget-title="Community"]` is absent from the DOM, the other widgets remain.
**Postconditions / cleanup:** restore Community via "Customize widgets..." → `ignored:false`.
> Auto-note: a synthetic `mouseenter` does not enable the icon's visibility — use Playwright `.hover()`; the close click itself works even without hover (fallback `click({force:true})`).

### Widgets-Controls-02 — A closed widget stays hidden after reload (persistence)
**Type:** regression
**Preconditions:** Home page open, Community visible.
**Steps:**
1. Close **Community** (close icon).
2. **Wait** for settings to be saved (async — see auto-note).
3. Reload the page (F5).
**Final check:**
- ✅ after reload **Community is absent**;
- ✅ `grok.userSettings.get('widgets')` → `communityWidget` = `{"ignored":true}`.
**Postconditions / cleanup:** restore Community via Customize → `ignored:false`, re-verify the reload.
> Auto-note: without a pause before the reload the case flakes (the `ignored:true` save does not reach the server). Captured live 2026-06-08.

---

## 4. Customize — "Customize widgets..." (Context Panel)

### Widgets-Customize-01 — Settings form: hiding and showing a widget
**Type:** functional
**Preconditions:** Home page open, Community visible.
**Steps:**
1. Click the **"Customize widgets..."** link under the widgets → reaction: a form of bool toggles opens in the Context Panel (one toggle per widget, label = friendlyName; visible ones are checked).
2. Uncheck **Community** → reaction: the widget disappears from the Home page.
3. Check **Community** again → reaction: the widget reappears.
**Final check:**
- ✅ the form has toggles for at least **Spotlight, Community, Usage, Reports**;
- ✅ after unchecking, `[widget-title="Community"]` is absent, `communityWidget` = `{"ignored":true}`;
- ✅ after re-checking, the widget is visible, `communityWidget` = `{"ignored":false}`.
**Postconditions:** state restored.
> Auto-note: the same toggle restores a widget previously closed via the close icon (`Widgets-Controls-01`) — a separate case is not needed, the mechanism is identical.

### Widgets-Customize-02 — Visibility setting persists after reload
**Type:** regression
**Preconditions:** Home page open.
**Steps:**
1. In Customize uncheck **Community**, wait for the save.
2. Reload the page.
**Final check:**
- ✅ after reload Community is hidden;
- ✅ the Customize form still shows Community unchecked.
**Postconditions / cleanup:** re-check the toggle → `ignored:false`.

---

## 5. Search — interaction between search and widgets

### Widgets-Search-01 — Entering a query hides the widgets, clearing brings them back
**Type:** functional
**Preconditions:** Home page open, the widgets panel visible.
**Steps:**
1. Type `aspirin` into the search box → (debounce ~500 ms) reaction: `.power-pack-widgets-panel` hides, `.power-pack-search-host` shows with results.
2. Clear the search box → reaction: the widgets panel is visible again, the results are hidden.
**Final check:**
- ✅ on input: the widgets panel is hidden, the results host is visible, path = `search?q=aspirin`;
- ✅ on clear: the widgets panel is visible, path = `search` (without `?q=`).
**Postconditions:** —

---

## 6. Navigation — navigating to the Home page

### Widgets-Nav-01 — The Home icon returns to Home, the widget set is stable
**Type:** smoke
**Preconditions:** A different view is open (e.g. a table / Browse node).
**Steps:**
1. Click the **Home** icon (`[name="icon-home"]`) → reaction: the Home page opens with the widgets.
2. Reload the page.
**Final check:**
- ✅ `.power-pack-welcome-view` with widgets is present;
- ✅ the set of visible widgets is identical before/after reload (Spotlight, Community, Usage, Reports).
**Postconditions:** —

---

## 7. Spotlight — multi-tab widget

### Widgets-Spotlight-01 — Rendering with tabs and "Demo of the day"
**Type:** smoke
**Preconditions:** Home page open, Spotlight visible.
**Steps:**
1. Inspect **Spotlight** → the tabs **Workspace, Spotlight, Favorites, Notifications, My Activity, Learn** and the **"💡 Demo of the day: <name>"** link at the bottom are visible.
**Final check:**
- ✅ all 6 tabs are present (`.d4-tab-header`);
- ✅ the **Notifications** tab shows an unread-count badge (if any);
- ✅ the "Demo of the day" link is present and leads to a demo when clicked.
**Postconditions:** —

### Widgets-Spotlight-02 — Switching tabs (incl. the Learn sub-tabs)
**Type:** functional
**Preconditions:** Spotlight visible.
**Steps:**
1. Click **Workspace** → pinned entities, the hint "Select a pinned item...".
2. Click **Spotlight** → Shared with me / Recent entities.
3. Click **Favorites** → favorites.
4. Click **My Activity** → the user's activity.
5. Click **Learn** → the sub-tabs **VIDEO / WIKI / DEMO / TUTORIALS** and topics (Meetings, Develop, Cheminformatics, Visualize, Explore); switch the sub-tabs.
**Final check:**
- ✅ each tab switches the `.power-pack-widget-content` content without errors;
- ✅ on the Learn tab all 4 sub-tabs are present, switching works.
**Postconditions:** —

### Widgets-Spotlight-03 — The Notifications tab and "Mark all as read"
**Type:** functional
**Preconditions:** Spotlight visible, there are unread notifications.
**Steps:**
1. Open the **Notifications** tab → reaction: a list of notifications, the **"Mark all as read"** link, an "N unread" counter.
2. Click **"Mark all as read"** → reaction: the unread counter resets to zero / the badge disappears.
**Final check:**
- ✅ notifications are displayed;
- ✅ after "Mark all as read" the unread badge is reset.
**Postconditions:** —
> Auto-note: the "there are unread" precondition is unstable — create a notification via a second user = share something or send a message, whichever is easier.

---

## 8. Community — links widget

### Widgets-Community-01 — Rendering the list of links
**Type:** smoke
**Preconditions:** Home page open, Community visible.
**Steps:**
1. Inspect the **Community** widget → a list of links (Platform Releases, Plugin releases, Macromolecules updates, Cheminformatics updates, Welcome to the Datagrok community, etc.).
**Final check:**
- ✅ there are at least a few links;
- ✅ each link points to `community.datagrok.ai` (check the `href`, not an actual navigation — external resource; but if it can be verified, that's good).
**Postconditions:** —

---

## 9. Usage and Reports — UsageAnalysis widgets (require permissions)

> Visible only to users in the **Developers / Administrators** groups (`canView`).
> The test user `opavlenko+playwright` has the administrator role → both widgets are visible.

### Widgets-Usage-01 — Rendering Usage and "Open Usage Analysis"
**Type:** functional
**Preconditions:** Home page under a user with permissions; Usage visible.
**Steps:**
1. Inspect the **Usage** widget → the **Users** and **Errors** charts (line charts on canvas), the **System** block with links (Datlas, Datlas DB, Grok Connect, Jupyter, Grok Spawner), the **"Open Usage Analysis"** link.
2. Click **"Open Usage Analysis"** → reaction: the UsageAnalysis application opens.
**Final check:**
- ✅ both charts, the System links, and "Open Usage Analysis" are present;
- ✅ the Usage Analysis application opens without errors.
**Postconditions / cleanup:** return to the Home page.

### Widgets-Reports-01 — Rendering Reports and "Open Reports"
**Type:** functional
**Preconditions:** Home page under a user with permissions; Reports visible.
**Steps:**
1. Inspect the **Reports** widget → the **"Open Reports"** link and a list of recent error reports.
2. Click **"Open Reports"** → reaction: a view with the list of reports opens.
**Final check:**
- ✅ the "Open Reports" link is present, the widget scrolls through recent reports (or is empty);
- ✅ the reports list opens without errors.
**Postconditions / cleanup:** return to the Home page.
> ⚠️ **Observation (2026-06-08):** on dev, Reports showed 4 entries "NullError: method not found: 'clientWidth' on null" —
> — recent auto-error-reports from the original error (likely the widget rendering in a zero-width container). The widget
> works (scrolls through reports), but the error itself can be ignored for now.

---

## 10. Permissions — gating by `canView`

### Widgets-Perm-01 — Usage/Reports hidden for a user without permissions
**Type:** negative
**Preconditions:** Home page under a user **without** the Developer/Administrator roles —
the secondary account `DATAGROK_SHARING_LOGIN` (`opavlenko+pwsharing`, storageState `.auth-sharing.json`).
**Steps:**
1. Under the non-admin account, open the Home page.
**Final check:**
- ✅ **Usage** and **Reports** are **not displayed**;
- ✅ only Spotlight and Community are visible.
**Postconditions:** —
> Auto-note: make sure the sharing account has no Developers/Administrators roles; otherwise set up a separate user.

### Widgets-Perm-02 — Usage/Reports visible to Developer/Administrator
**Type:** functional
**Preconditions:** Home page under the main account `opavlenko+playwright` (has the admin role).
**Steps:**
1. Open the Home page.
**Final check:**
- ✅ all 4 are present in DOM order: Spotlight, Community, Usage, Reports.
**Postconditions:** —

---

## 11. Tag summary (for Playwright)

- `@smoke`: Display-01, Nav-01, Spotlight-01, Community-01.
- `@functional`: Controls-01, Customize-01, Search-01, Spotlight-02/03, Usage-01, Reports-01, Perm-02.
- `@regression`: Controls-02, Customize-02.
- `@negative`: Perm-01.

**Total: 15 cases.**
