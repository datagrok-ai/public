# Dashboards and sharing

Playlist: **Dashboards and sharing**

## Plan

What the platform does today (core `project_meta.dart` `ProjectMeta.publish`, `shell_project.dart` `openProject`,
help `datagrok/concepts/project/dashboard.md`):

- A dashboard is a project: tables plus their views and layout. The ribbon's **SAVE** opens the **Save project**
  dialog (name, description, picture, "Save to a space (optional)", Presentation mode, the list of tables with a
  per-table **Data sync** switch and **Creation script**). A new dashboard is private: right after the first save
  the **Share** dialog opens, and nobody else sees it until it is shared. New dashboards go to the user's own
  space: Browse > **My stuff** (My dashboards, Recent) and Browse > **Dashboards** (a gallery of every dashboard
  the user can open, searchable by name or #tag).
- Saving a dashboard that is already on the server offers three choices (radio buttons at the top of the dialog):
  - **Save original project** — overwrites it for everyone ("Save all data and visualizations"). Unavailable
    without the Edit privilege.
  - **Save a copy** — a new project ("Copy of ..."); for each table, Move, **Link** (keeps showing the original's
    data) or **Clone** (an independent copy). The original is not changed.
  - **Save personal view customizations** — the platform's name for changing a dashboard for yourself only: the
    views and their arrangement go into the user's own settings (`projectCustomViews`/`projectCustomLayout`),
    the name is greyed out, the project on the server is not touched. On the next open, a balloon says "Project
    has personal view customizations"; the project's **Custom views** pane in the context panel says so too and
    has a **Reset** button that goes back to the shared version.
- Sharing: right-click a dashboard > **Share...** (or the context panel's Sharing pane): users, groups or emails,
  an access level (View and use, Edit, Full access...), **Send notifications**. The recipient finds it under My stuff
  > Shared with me > <owner>. The URL (Links...) opens it only for someone who already has access.
- Layouts: Toolbox > **Layouts** > Save stores just the arrangement of viewers; its card is offered for any table
  with matching columns, and a click applies it.
- Data: per table, **Data sync** on = the creation script runs again on open (fresh data from the source); off =
  the data is stored with the dashboard (a snapshot).

Clips, most important first (one idea each):

1. Save a dashboard and find it again (Save, the Share offer, Browse > Dashboards, search, My stuff > Recent,
   Add To Favorites, open from the Favorites panel).
2. Save changes to the original or as a copy (Save original project, Save a copy, the original unchanged).
3. Personal view customizations (change a dashboard for yourself, the balloon on reopen, Custom views > Reset).
4. Share a dashboard (Share dialog, access levels, notifications, the link; the other user finds it in their
   Browse > Dashboards and can only save a copy. My stuff > Shared with me was left out: on the shared stand it
   lists other features' leftovers).
5. Layouts (save the arrangement in Toolbox > Layouts, apply it to another table).
6. Where a dashboard's data comes from (Data sync on/off, the creation script).

Added after the six:

7. Share a space with your team (Move to Space..., share the space once, access inherited by what is inside).
8. Keep a dashboard in Favorites: merged into 1 (a 29 s clip was too short on its own).
9. Presentation mode (planned, not delivered): the Save dialog's Presentation mode switch is not kept on this
   build, so the dashboard reopens in design mode. The feature is written
   (`08-presentation-mode.feature`) and fails on that step until the product is fixed.

Everything a clip saves (projects, layouts, shares, personal customizations) is deleted at the end of its feature
and swept at its start.

## 01 - Save a dashboard and find it again (1:12)

Title: Save a dashboard in Datagrok and find it again

Description: A table with a scatter plot, a histogram and a bar chart is saved as a dashboard with the SAVE button
at the top of the table view. Right after the first save Datagrok offers to share it; until you do, only you can
see it. After closing it, you find it in Browse > Dashboards by searching for its name, see who made it in the
Context Panel, and find it under My stuff > Recent too. Right-click it and pick Add To Favorites (only for you, or
for a group), and the star on the left sidebar lists it, ready to open with a double-click and all its viewers.
Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/06-dashboards-and-sharing/01-save-and-find-a-dashboard.feature

## 02 - Save changes to the original or as a copy (1:12)

Title: Save a dashboard to the original or as a copy in Datagrok

Description: When you save a dashboard that is already on the server, the Save project dialog asks where the
change goes. "Save original project" updates the dashboard for everyone who opens it; "Save a copy" makes a new
dashboard under its own name and leaves the original as it was. The video does both, closes the dashboard, and
opens the original again from Browse > Dashboards to show which change it kept. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/06-dashboards-and-sharing/02-save-to-the-original-or-as-a-copy.feature

## 03 - Change a dashboard just for yourself (1:12)

Title: Change a shared dashboard just for yourself in Datagrok

Description: To change a dashboard only for yourself, pick "Save personal view customizations" in the Save project
dialog. Your version is kept in your own settings, the dashboard on the server stays as it was, and when you open
it again from Browse > Dashboards Datagrok shows your version. To go back, click the dashboard in Browse >
Dashboards and use Reset in the Custom views section of the Context Panel; opened again, it is the shared version.
Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/06-dashboards-and-sharing/03-personal-view-customizations.feature

## 04 - Share a dashboard (1:11)

Title: How to share a dashboard in Datagrok

Description: A new dashboard is visible only to you until you share it. Right-click it in Browse > Dashboards and
pick Share..., then add a user, a group or an email and choose the access level (View and use, Edit, Full access).
In the Context Panel, Links... gives the dashboard's address, which opens only for people who have access, and the
Sharing section shows who they are. Then, signed in as the colleague, the video finds the dashboard in their own
Browse > Dashboards; with View and use they can open it and save a copy, but not save over the original. Uses the
demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/06-dashboards-and-sharing/04-share-a-dashboard.feature

## 05 - Save a layout and apply it to another table (0:46)

Title: Save a layout of viewers and reuse it in Datagrok

Description: A layout keeps the viewers and how they are set up, without the data. Save it from the Toolbox, in the
Layouts section, and rename it in the Context Panel. Any other table with the same columns, here demog-1000.csv
opened from Browse > Files > Demo, then lists the layout in its own Layouts section, and one click arranges that
table the same way. Uses the demog and demog-1000 demo datasets.

Feature: public/packages/Tutorials/bdd/features/howto/06-dashboards-and-sharing/05-save-and-apply-a-layout.feature

## 06 - Where a dashboard gets its data (0:42)

Title: Dashboard data in Datagrok: kept in sync or saved as a snapshot

Description: When you save a dashboard, the Save project dialog lists each table with a Data sync switch. With Data
sync on, the dashboard runs the table's creation script every time it opens, so it reads the file or query again
and shows current data; with it off, the data is stored with the dashboard as it is now. The video opens demog.csv
from Browse > Files > Demo, looks at its creation script, and saves the dashboard as a snapshot.

Feature: public/packages/Tutorials/bdd/features/howto/06-dashboards-and-sharing/06-where-the-data-comes-from.feature

## 07 - Share a space with your team (1:00)

Title: Share a space with your team in Datagrok

Description: Instead of sharing dashboards one at a time, put them in a space and share the space. The video moves
a dashboard into a space with Move to Space... (right-click it in Browse > Dashboards), expands Browse > Spaces and
opens the space, and shares the space with a colleague. The dashboard's own Share dialog then shows the same
access, inherited from the space. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/06-dashboards-and-sharing/07-share-a-space.feature