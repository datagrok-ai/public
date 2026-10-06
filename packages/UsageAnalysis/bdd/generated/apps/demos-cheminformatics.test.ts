/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/demos-cheminformatics.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {customFiredWith, listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, packageInstalled, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {moreHighlight, noBalloons, noErrors, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Cheminformatics demos open from Browse > Apps > Demo with their content", () => {
  const session = feature(test, "features/apps/demos-cheminformatics.feature", import.meta.url);
  test("The SAR Matrix demo opens with its SAR Matrix Viewer viewer and its table [demo=SAR Matrix, node=SAR-Matrix, viewer=SAR Matrix Viewer, rows=10000]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(23, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(24, "When user clicks on Apps---Demo---Cheminformatics---SAR-Matrix tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---SAR-Matrix tree node inside browse tree")));
    await session.step(25, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | SAR Matrix\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | SAR Matrix"));
    await session.step(26, "And the \"SAR Matrix\" view should be current", () => viewIsCurrent(page, "SAR Matrix"));
    await session.step(27, "And SAR Matrix Viewer viewer should be visible", () => shouldBe(page, el("SAR Matrix Viewer viewer"), "visible"));
    await session.step(28, "And the table should have 10000 rows", () => rowCount(page, 10000));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Chemical Space demo opens with its scatter plot viewer and its table [demo=Chemical Space, node=Chemical-Space, viewer=scatter plot, rows=1000]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(23, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(24, "When user clicks on Apps---Demo---Cheminformatics---Chemical-Space tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---Chemical-Space tree node inside browse tree")));
    await session.step(25, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | Chemical Space\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | Chemical Space"));
    await session.step(26, "And the \"Chemical Space\" view should be current", () => viewIsCurrent(page, "Chemical Space"));
    await session.step(27, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(28, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Molecule Activity Cliffs demo opens with its scatter plot viewer and its table [demo=Molecule Activity Cliffs, node=Molecule-Activity-Cliffs, viewer=scatter plot, rows=200]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(23, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(24, "When user clicks on Apps---Demo---Cheminformatics---Molecule-Activity-Cliffs tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---Molecule-Activity-Cliffs tree node inside browse tree")));
    await session.step(25, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | Molecule Activity Cliffs\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | Molecule Activity Cliffs"));
    await session.step(26, "And the \"Molecule Activity Cliffs\" view should be current", () => viewIsCurrent(page, "Molecule Activity Cliffs"));
    await session.step(27, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(28, "And the table should have 200 rows", () => rowCount(page, 200));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The R-Group Analysis demo opens with its trellis plot viewer and its table [demo=R-Group Analysis, node=R-Group-Analysis, viewer=trellis plot, rows=200]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(23, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(24, "When user clicks on Apps---Demo---Cheminformatics---R-Group-Analysis tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---R-Group-Analysis tree node inside browse tree")));
    await session.step(25, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | R-Group Analysis\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | R-Group Analysis"));
    await session.step(26, "And the \"R-Group Analysis\" view should be current", () => viewIsCurrent(page, "R-Group Analysis"));
    await session.step(27, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
    await session.step(28, "And the table should have 200 rows", () => rowCount(page, 200));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Matched Molecular Pairs demo opens with its Matched Molecular Pairs Analysis viewer and its table [demo=Matched Molecular Pairs, node=Matched-Molecular-Pairs, viewer=Matched Molecular Pairs Analysis, rows=20267]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(23, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(24, "When user clicks on Apps---Demo---Cheminformatics---Matched-Molecular-Pairs tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---Matched-Molecular-Pairs tree node inside browse tree")));
    await session.step(25, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | Matched Molecular Pairs\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | Matched Molecular Pairs"));
    await session.step(26, "And the \"Matched Molecular Pairs\" view should be current", () => viewIsCurrent(page, "Matched Molecular Pairs"));
    await session.step(27, "And Matched Molecular Pairs Analysis viewer should be visible", () => shouldBe(page, el("Matched Molecular Pairs Analysis viewer"), "visible"));
    await session.step(28, "And the table should have 20267 rows", () => rowCount(page, 20267));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Similarity & Diversity Search demo opens with its Chem Similarity Search viewer and its table [demo=Similarity & Diversity Search, node=Similarity-&-Diversity-Search, viewer=Chem Similarity Search, rows=1000]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(23, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(24, "When user clicks on Apps---Demo---Cheminformatics---Similarity-&-Diversity-Search tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---Similarity-&-Diversity-Search tree node inside browse tree")));
    await session.step(25, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | Similarity & Diversity Search\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | Similarity & Diversity Search"));
    await session.step(26, "And the \"Similarity & Diversity Search\" view should be current", () => viewIsCurrent(page, "Similarity & Diversity Search"));
    await session.step(27, "And Chem Similarity Search viewer should be visible", () => shouldBe(page, el("Chem Similarity Search viewer"), "visible"));
    await session.step(28, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Scaffold Tree demo opens with its Scaffold Tree viewer and its table [demo=Scaffold Tree, node=Scaffold-Tree, viewer=Scaffold Tree, rows=1000]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(23, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(24, "When user clicks on Apps---Demo---Cheminformatics---Scaffold-Tree tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---Scaffold-Tree tree node inside browse tree")));
    await session.step(25, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | Scaffold Tree\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | Scaffold Tree"));
    await session.step(26, "And the \"Scaffold Tree\" view should be current", () => viewIsCurrent(page, "Scaffold Tree"));
    await session.step(27, "And Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
    await session.step(28, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Rows selected in the Chemical Space demo light up in its scatter plot", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(20, "And Apps---Demo---Cheminformatics tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Cheminformatics tree node inside browse tree")));
    await session.step(44, "Given user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(45, "When user clicks on Apps---Demo---Cheminformatics---Chemical-Space tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Cheminformatics---Chemical-Space tree node inside browse tree")));
    await session.step(46, "Then the \"demo-loaded\" custom event should have fired with path \"Cheminformatics | Chemical Space\"", () => customFiredWith(page, "demo-loaded", "path", "Cheminformatics | Chemical Space"));
    await session.step(47, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
    await session.step(48, "And user selects the first 100 rows", () => selectFirstRows(page, 100));
    await session.step(49, "Then 100 rows should be selected", () => selectedRowCount(page, 100));
    await session.step(50, "And scatter plot viewer should show more selection highlight than before", () => moreHighlight(page, el("scatter plot viewer")));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
});
