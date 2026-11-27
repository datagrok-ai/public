# Review: 08-end-to-end (From data to answers)

Fresh-eyes review (2026-11-22, ~21:45-22:00): feature lines, descriptions.md, and contact sheets at the key lines and
results of each clip (12-16 frames for these longer clips). Every delivered video matched its latest guide take.
The m:ss in descriptions.md are the delivered lengths.

| NN | Clip | Verdict |
|---|---|---|
| 01 | From a file to a shared dashboard | OK. Every move is filmed, from Browse > Files > Demo to the Dashboards gallery. Needs the lead: (1) the unfiltered Dashboards gallery shows stand leftovers for ~3 s before the search narrows it ("Auto Test- Formula Refresh", "Auto Test- Add New Col", "cicd_tests", several *_demo projects); (2) the colleague the dashboard is shared with is the fixture user "bddsecond", typed and shown on screen: a test-looking name (the same in the 06 section). |
| 02 | Explore a new dataset | refilmed (1:12): the on-screen claim "Horsepower goes with price most" was false (engine size r 0.91, horsepower and curb weight 0.83), "the most powerful" was false for the Mercedes (123 hp; mercury 175), and with all 17 columns the correlation plot ran off the right edge of its dock, so the cells the lines talk about sat at the cut edge. Now cars is cut to 8 columns off camera at the start (status bar line "30 rows, 8 columns"), the plot is 7x7 and whole, the lines point at engine size, and the scatter plot is engine size against price (the three selected cars have the three biggest engines, at the top right). "What makes them dear" became "Why so expensive". Title "...in a minute" became "Datagrok: a first look at a new dataset". The view's top-right label reads "null" after the double-click (known product bug of the correlation plot). |
| 03 | Compare groups | text fixed (re-rendered from the same take, 1:16): "The p-value says it is not chance" overclaimed; now "Men weigh more; the tiny p-value says it is unlikely to be chance". Title was a question ("do two groups differ?"): now "Datagrok: compare two groups with a filter, a box plot and a layout"; the description's opening question is a statement now. The Layouts pane shows two stand layouts named "Table" (leftovers). |
| 04 | Clean and reshape a table | OK. |
| 05 | Combine two tables to answer a question | text fixed: the description opened with a question ("Which products sell the most?"); now a statement. The answer (Camembert Pierrot, 1,577 units) matches the Northwind data. Result table named "result" (platform default). |
| 06 | From molecules to a shortlist | OK. The ribbon shows Bio next to Chem on spgi (known product bug). |

Verdicts: 3 OK, 2 text fixed (03 re-rendered), 1 refilmed (02), 0 need work; 01 waits on the lead's stand cleanup
and sharing-user decision.

## Titles and descriptions

Each clip states its task and walks through it in plain sentences; the series' cross-links to the single-idea
playlists are useful. Two questions (03 title, 05 opening) and one hype phrase ("in a minute", 02) removed.

## Jumps

None after the start. The film.sh warning scan is clean except 02's new column removals, which are part of the
starting state (before the first visible step). Views, panels, viewers and dialogs all come from filmed clicks;
Browse is expanded on camera in 01.

## Systematic

- End-to-end clips state facts about the data on screen ("goes with price most", "the most powerful", "not chance");
  these deserve a check against the numbers in the frame, as the clip itself is the only proof the viewer gets.
  02 had two false ones. Worth a line in the brief: verify every data claim against the values on screen.
- A viewer that compares every numeric column (correlation plot) needs a narrow table, or it runs off its dock at
  1080p with the toolbox and a second viewer open.
- Galleries (Dashboards, Layouts) show stand leftovers whenever a clip opens them unfiltered.
