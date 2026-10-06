---
feature: projects
realizes_atlas: []
realizes: []
priority: p2
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  The local-machine cases go through the OS file picker (not automatable);
  the remaining cases are gaps in the automated uploading-spec.ts source
  matrix, kept manual until folded into the spec.
automation_candidate: true
related_bugs: []
---

# Project uploading — manual companion

Cases that are absent from the automated matrix in `uploading.md`. Follow the
same save/verify/cleanup pattern as `uploading.md`: save, close all, reopen
from Browse, verify content, then delete the project.

## Project from two local files

1. Open two files from the local machine (drag onto the platform or
   File > Open)
2. Save them as one project named `UI_Local_NoSync`
3. Reopen the project — `UI_Local_NoSync` opens the stored snapshot tables
4. Delete the created project

## Project from an SDF file

1. Open `System:AppData/Chem/mol1K.sdf` from Files
2. Save it as a project named `UI_Sdf_Sync` with Data Sync on
3. Reopen the project — molecules render, no errors
4. Delete the created project

## Project from the Get Top 100 query result

1. In Browse > Databases, run **Get Top 100** on a demo DB table (e.g.
   Northwind `orders`)
2. Save the result as a project — `UI_Top100_Sync` with Data Sync on, then
   `UI_Top100_NoSync` with it off
3. Reopen each project — with Sync on the query re-executes; with Sync off
   the stored result opens
4. Delete both created projects

## Upload dialog from the Scratchpad

1. Open any table; open the Dashboards panel (from the Sidebar) and click the **Save** button near the New Dashboard project
2. In the dialog: enter the name `UI_Scratchpad_Upload`, enter a description,
   turn on the **Presentation mode** switch, click OK
3. In the Share dialog that opens, click OK
4. Reopen the project — the name and description are shown on its card, and
   opening it starts in presentation mode
5. Delete the created project

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/mol1K.sdf"]
}
