/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/files/newick-file.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [dendrogram.cp.newick-file-open-via-files-browser]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openNewickFile, previewTreeLeaves, tableOfTreeLeaves} from '../../bindings/tree-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, fillsParent, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableTagIsFile} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A newick file opened from Browse and through its file handler", () => {
  const session = feature(test, "features/files/newick-file.feature", import.meta.url);
  test("A newick file opened from Browse and through its file handler", {tag: ["@journey", "@realizes:dendrogram.cp.newick-file-open-via-files-browser"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(20, "And \"Files > App Data\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data\" tree node inside browse tree")));
    await session.step(21, "And \"Files > App Data > Dendrogram\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data > Dendrogram\" tree node inside browse tree")));
    await run.scenario("The file handler makes a tree table with a Dendrogram of the same leaves", async () => {
      await session.step(24, "When user opens the newick file \"System:AppData/Dendrogram/data/nwk1.nwk\" with its file handler", () => openNewickFile(page, "System:AppData/Dendrogram/data/nwk1.nwk"));
      await session.step(25, "Then Dendrogram viewer should be visible", () => shouldBe(page, el("Dendrogram viewer"), "visible"));
      await session.step(26, "And the table should have tag \".newick\" equal to the text of \"System:AppData/Dendrogram/data/nwk1.nwk\" file", () => tableTagIsFile(page, ".newick", "System:AppData/Dendrogram/data/nwk1.nwk"));
      await session.step(27, "And the \"leaves\" reading of Dendrogram viewer should be \"leaf1, leaf2, leaf3\"", () => readingReads(page, "leaves", el("Dendrogram viewer"), "leaf1, leaf2, leaf3"));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on the file previews its tree whether or not PhyloTreeViewer is installed", async () => {
      await session.step(32, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(33, "When user clicks on \"Files > App Data > Dendrogram > data\" tree node inside browse tree", () => clickOn(page, el("\"Files > App Data > Dendrogram > data\" tree node inside browse tree")));
      await session.step(34, "And user clicks on nwk1.nwk link in gallery", () => clickOn(page, el("nwk1.nwk link in gallery")));
      await session.step(35, "Then the newick preview should hold a tree with leaves \"leaf1, leaf2, leaf3\"", () => previewTreeLeaves(page, "leaf1, leaf2, leaf3"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click on the file shows its tree", async () => {
      await session.step(39, "When user clicks on \"Files > App Data > Dendrogram > data\" tree node inside browse tree", () => clickOn(page, el("\"Files > App Data > Dendrogram > data\" tree node inside browse tree")));
      await session.step(40, "Given the \"PhyloTreeViewer\" package is installed", () => packageInstalled(page, "PhyloTreeViewer"));
      await session.step(41, "When user double-clicks on nwk1.nwk link in gallery", () => doubleClickOn(page, el("nwk1.nwk link in gallery")));
      await session.step(42, "Then PhylocanvasGL viewer should be visible", () => shouldBe(page, el("PhylocanvasGL viewer"), "visible"));
      await session.step(43, "And PhylocanvasGL viewer should fill its parent", () => fillsParent(page, el("PhylocanvasGL viewer")));
      await session.step(44, "And the table of PhylocanvasGL viewer should hold a tree with leaves \"leaf1, leaf2, leaf3\"", () => tableOfTreeLeaves(page, el("PhylocanvasGL viewer"), "leaf1, leaf2, leaf3"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on the file previews the same tree", async () => {
      await session.step(48, "When user clicks on \"Files > App Data > Dendrogram > data\" tree node inside browse tree", () => clickOn(page, el("\"Files > App Data > Dendrogram > data\" tree node inside browse tree")));
      await session.step(49, "And user clicks on nwk1.nwk link in gallery", () => clickOn(page, el("nwk1.nwk link in gallery")));
      await session.step(50, "Then PhylocanvasGL viewer should be visible", () => shouldBe(page, el("PhylocanvasGL viewer"), "visible"));
      await session.step(51, "And PhylocanvasGL viewer should fill its parent", () => fillsParent(page, el("PhylocanvasGL viewer")));
      await session.step(52, "And the table of PhylocanvasGL viewer should hold a tree with leaves \"leaf1, leaf2, leaf3\"", () => tableOfTreeLeaves(page, el("PhylocanvasGL viewer"), "leaf1, leaf2, leaf3"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
