### The checklist for testing Spaces

1. **Create**  
   - Root  
   - Child  
   - Multiple levels of nesting  

2. **Update**  
   - Rename root  
   - Rename child  
   - Add to favorites / Remove from favorites  
   - Copy / Rename copy  
   - Link / Rename link  
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
   - Linked files, projects, and functions from a space  
   - Linked space  
   - Space (root and child)  
   - Files in a space  
   - Projects in a space  

7. **Forbidden**  
   - Moving functions from packages to spaces  
   - Dragging a root space into its child (including linked spaces located within the child)  
   - Dragging a child space to the root *Spaces* level (waiting for clarification)  
   - Creating duplicate space names on the same level  
   - Using an empty name (`""`)
---
{
  "order": 4
}