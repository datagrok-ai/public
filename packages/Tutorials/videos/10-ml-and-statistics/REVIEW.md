# Review: 10-ml-and-statistics (Machine learning and statistics)

Fresh-eyes review (2026-11-22, ~21:45-22:15) of the 9 delivered clips: feature, descriptions.md block, and contact
sheets at the key lines and results of each delivered take (every delivered video matched its latest take's length
before the review). Judged on: is the one idea clear to someone new to Datagrok, is each line true when it shows,
plain wording, jumps. The m:ss in descriptions.md match the delivered videos.

| NN | Clip | Verdict |
|---|---|---|
| 01 | What the ML menu does | text fixed: the description said "the later videos of this playlist show each tool", but Impute, k-means, DBSCAN, MCL and the t-test have no video here; now "most of these tools have a video of their own in this playlist". A menu tooltip lags one item behind for a moment (Impute's tooltip beside Pareto Front at 28 s); transient. |
| 02 | Reduce many columns with PCA | OK. Red and white wines sit apart along PC1 in the take, as the line says. |
| 03 | Map similar rows with UMAP | OK, text added. The maker's question, checked in the delivered take: red wines fill the lower tail of the map, white ones the rest, a few strays each way, so "land in different parts of the map" is true. The description now says the map's shape changes from run to run while the neighbours don't. |
| 04 | Train a model | captions fixed, re-rendered from the same take (1:05). "Confusions: which species it mixes up" sat over a perfectly diagonal 50/50/50 matrix: Forge's Confusions panel shows the model trained on all rows predicting those same rows (validation accuracy 0.927 is right beside it), so a newcomer reads "it mixes up nothing". Now "Confusions: predicted against actual, on the rows it learned from". "Method: XGBoost here; every method runs in the browser" overclaimed (Forge has server methods, e.g. Chemprop for molecules) and is now "training runs in your browser". The description says the matrix is on the training rows. Models sweep untouched (no refilm). |
| 05 | Apply a model to new data | OK, text added: the description now says the model is the one from the previous video and names the iris dataset. The line "models that fit this table are listed first" is true but the list is never opened on screen (one model on the stand); left as is. |
| 06 | Find what drives a value with multivariate analysis | refilmed (1:12). The maker's question: the Loadings and Scores lines sat at the top over the Observed vs. Predicted title they don't describe; now `the video says ... at the bottom`. Added "Leave out price itself: it is the value to explain" on the search, since price was unticked from Using without a word. "Explained Variance: how much of price each component explains" misread the chart (the bars are cumulative: 1, 2, 3 components explain 0.70, 0.82, 0.92 of price); now "share of price explained by 1, 2 and 3 components". The bottom Loadings line clips the "horsepower" label at the chart's lower right for ~2.5 s. |
| 07 | Compare groups with ANOVA | text fixed: the title "One-way ANOVA in Datagrok: do groups differ" read as a question; now "compare a value across groups". The description says the test is about the average. |
| 08 | Find the best trade-offs with a Pareto front | OK. Checked against data/demo/cars.csv: the 11 green cars are exactly the price/horsepower Pareto front (toyota1, mitsubishi, chevrolet, isuzu, subaru, dodge1, toyota3, alfaromeo, mercury, jaguar, porsche). |
| 09 | Cluster rows into a tree | refilmed (1:05). The maker's question: the Cluster column was colored by an off-camera Given right after Assign, so it seemed to come out of Assign already painted (in 6 px rows an uncolored column looks empty). That was a jump; now it is a filmed Color Coding > Categorical from the column header's menu, then the menu is closed (the pick is a radio item and leaves it open). Species is colored off camera with the same default colors in the clusters' order (setosa blue, virginica light orange, versicolor green), so matching colors mean agreement. New closing line "Setosa gets a cluster of its own; the other two species overlap", true in the take: cluster 1 is all 50 setosa, cluster 2 virginica only, cluster 3 versicolor plus the remaining virginica. The 6 px rows (starting state) leave the numeric columns, Medoid Rank and Avg Distance looking empty; the description explains the short rows. |

Verdicts: 3 OK (02, 05 with a description line added, 08), 2 text fixed (01, 07), 1 OK with a sentence added (03),
1 captions fixed by re-render (04), 2 refilmed (06, 09).

## Data claims (every on-screen line and description sentence that states a result, checked against the frame)

All true after the fixes: PCA's red/white split along PC1 (02); the UMAP map's red tail (03); Performance's 5-fold
validation and the two Runs (04); all 9 predicted species match the true species of iris rows 1, 51, 101, 11, 61,
111, 21, 81, 121 (05); fifteen specs, Using (15), eng.size 1.52 and horsepower 1.37 on top of Variable Importance,
0.70 / 0.82 / 0.92 of price explained (06); six DIS_POP groups, F 44.44 above F-critical 2.26, p 7.31E-31 below 0.05
(07); the 11 green cars = the computed Pareto front (08); setosa alone in cluster 1 (09). The two misreadings
found were 04's Confusions line and 06's Explained Variance line, both fixed.

## Titles and descriptions

Plain and specific, no hype words; each names its menu path and demo dataset. Fixed: one question-shaped title (07),
one claim the playlist does not keep (01), one detail that made a panel look wrong (04).

## Systematic

- Starting states are set off camera (rows 6 px, column colors, the toolbox hidden); after the first visible step the
  only Givens left are `user watches the task bar` (03, 09), which changes nothing on screen (film.sh still prints
  its no-jumps warning for them).
- 04 was re-rendered, not refilmed: only two caption texts changed, so the feature and the take's steps.json got the
  same new lines and guide-render.py redrew the same take (the Forge model/run sweep was not rerun).
- A grid header's Color Coding picks are radio items that leave the context menu open; a feature that picks one
  needs `user closes the context menu` after it.
- Lines placed by the renderer can land over a neighbouring chart's title when a view is tiled with charts (06);
  a `the video says ... at the bottom` pin is the fix, since `# caption:` cannot be pinned.
- Product/design note for the lead (Forge): Confusions and Predicted vs Actual show the full-data model on its own
  training rows while Performance and Per-class show 5-fold validation; with XGBoost or Random Forest on iris the
  matrix is always a perfect diagonal next to a 0.927 validation accuracy. Out-of-fold predictions there would
  match the scores beside them.

## Merge round (2026-11-22, ~23:15-23:35)

The lead: clips under ~35 s are too thin. One merge and one extension, each filmed, reviewed by frames and delivered:

| NN | Clip | What changed |
|---|---|---|
| 04 | Train a model and apply it (1:34) | 04 "Train a model" (1:05) and 05 "Apply a model to new data" (0:33) in one take: everything 04 showed, then the "new flowers" view tab is clicked (both tables opened off camera at the start, "new flowers" first) and Apply Model picks the model just saved (claimed by name), matches the columns, adds Species (predicted). All nine predictions on screen match the true species of iris rows 1, 51, 101, 11, 61, 111, 21, 81, 121. The Apply fixture model (Forge:trainModel off camera) is gone: the clip's own saved Random Forest model is applied. One caption re-rendered from the take ("Model: the one just saved is already picked"). The model, its application and the runs were swept at start and end. |
| 08 | Find the best trade-offs with a Pareto front (0:48) | Extended from 0:32 with what the front means, on three cars, each pointed at with its tooltip in view: grey Mercedes (123 hp, 31,600), the Alfa Romeo that beats it (154 hp, 16,500), green Porsche (the dearest, 37,028, and the most powerful, 207 hp); then the PC plot. Checked against data/demo/cars.csv and claimed in the feature (hovered row, model, Pareto optimality, the maxima). One caption corrected by re-render: "every car across all its columns" was false (the PC plot shows 10 of the 15 numeric columns); now "one line per car across its columns, green on the front". |

Absorbed: 05-apply-a-model.feature (and its generated spec), "05 - Apply a model to new data.mp4"; superseded by the
new name: "04 - Train a model.mp4".
