# Initial UX — Implementation Status

## Spotlight widget

| Tab            | Status | Notes                                                         |
|----------------|--------|---------------------------------------------------------------|
| Spotlight      | Done   | Action required, Shared with me, Recent sections             |
| Workspace      | Done   | See details below                                            |
| Favorites      | Done   |                                                               |
| Notifications  | Done   | Unread count badge                                           |
| My Activity    | Done   | Deduplication of consecutive identical events                |
| Learn          | Done   |                                                               |

## Workspace tab

### Left pane — pinned entity list

| Feature                                   | Status | Notes                                           |
|-------------------------------------------|--------|-------------------------------------------------|
| "Myself only" section                     | Done   |                                                 |
| Per-group sections                        | Done   | One section per group the user belongs to       |
| "Continue where you left off" (recent)   | Done   | Up to 5 recent items with timestamps            |
| Pencil icon → Manage Favorites dialog     | Done   | Admin-only; full search + add/remove            |

### Pinning entities

| Method                                                | Status   | Notes                                          |
|-------------------------------------------------------|----------|------------------------------------------------|
| Drag-and-drop to Workspace widget                     | Done     | Drop zones per group + "Myself only"           |
| Right-click → "Group favorites" submenu               | Not done | No context menu integration yet                |

### Right pane — entity controls

| Entity type    | Controls shown                        | Status   |
|----------------|---------------------------------------|----------|
| Query          | Parameter editor + Run button         | Done     |
| Function       | Parameter editor + Run button         | Done     |
| Project/Dashboard | Project card (via ObjectHandler)   | Done     |
| App            | App header                            | Not done — shows generic "Open" hint |
| File / Share   | Details (type, path, modified date)   | Done                                 |
| DB Connection  | List of related queries               | Not done — shows generic "Open" hint |

### Bottom preview area

| Entity type    | Preview shown                         | Status   |
|----------------|---------------------------------------|----------|
| Query          | Results as grid (dynamically updated) | Done     |
| Function       | Results (DataFrame or custom view)    | Done     |
| Project        | Dashboard preview card                | Done     |
| App            | App preview                           | Not done |
| File / Share   | File/folder preview via ObjectHandler | Done     |
| DB Connection  | (none planned)                        | N/A      |

## Key source files

| File                                             | Purpose                                  |
|--------------------------------------------------|------------------------------------------|
| `src/spotlight/spotlight-widget.ts`              | Main widget, all tabs                    |
| `src/spotlight/workspace-tab.ts`                 | Workspace tab (list + controls + preview)|
| `src/spotlight/group-favorites.ts`               | Favorites API wrappers                   |
| `src/spotlight/manage-favorites-dialog.ts`       | Search & pin/unpin dialog                |
| `src/spotlight/preview-host.ts`                  | Preview area below the Spotlight widget  |
| `src/welcome-view.ts`                            | Integration with Welcome view            |
| `css/power-pack.css` (lines 1083–1630)           | All workspace/spotlight styles           |
