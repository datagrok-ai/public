### Spaces: scenarios not covered by server autotests

2. **Update**
   - Copy / Rename copy
   - Rename link
   - Rename entity in space (from both Browse Tree and view)

3. **Drag and Drop**
   - From local storage (group and single file) to a space (in both Browse Tree and view)
   - From one space to another (in both Browse Tree and view):
     - Link
     - Copy
     - Move
   - A group of files from space view to another space:
     - Link
     - Copy
     - Move
   * From Dashboards:
     - A project containing linked tables with viewers (with both **link** and **move** options) -  Verify that the linking between the tables is preserved after the move
   * From Files:
     * A **file/group of files** using each of the following options:
       - Link
       - Copy
       - Move
     * A **folder** using each of the following options:
       - Link
       - Copy
       - Move
   - From Packages
   - Alt + dragging — default action is link

4. **Read**
   - Preview entities (file, project, child space, linked file, linked project, linked function)
   - Open file or project
   - Search in the view

5. **Share**
   - Linked entity
   - Entity
   - Root space
   - Child space

6. **Delete** (in both Browse Tree and view)
   - Linked space

7. **Forbidden**
   - Moving functions from packages to spaces
   - Dragging a root space into its child (including linked spaces located within the child)
   - Dragging a child space to the root *Spaces* level (waiting for clarification)
---
{
  "order": 5
}