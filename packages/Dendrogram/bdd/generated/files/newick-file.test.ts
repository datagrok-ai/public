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
import {openNewickFile, tableOfTreeLeaves} from '../../bindings/tree-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, fillsParent, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableColumns, tableRows, tableTagIsFile} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A newick file opened from Browse and through its file handler", () => {
  const session = feature(test, "features/files/newick-file.feature", import.meta.url);
  test("A newick file opened from Browse and through its file handler", {tag: ["@journey", "@realizes:dendrogram.cp.newick-file-open-via-files-browser"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(13, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(14, "And \"Files > App Data\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data\" tree node inside browse tree")));
    await session.step(15, "And \"Files > App Data > Dendrogram\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data > Dendrogram\" tree node inside browse tree")));
    await run.scenario("A double-click on the file shows its tree", async () => {
      await session.step(18, "When user clicks on \"Files > App Data > Dendrogram > data\" tree node inside browse tree", () => clickOn(page, el("\"Files > App Data > Dendrogram > data\" tree node inside browse tree")));
      await session.step(19, "And user double-clicks on nwk1.nwk link in gallery", () => doubleClickOn(page, el("nwk1.nwk link in gallery")));
      await session.step(20, "Then PhylocanvasGL viewer should be visible", () => shouldBe(page, el("PhylocanvasGL viewer"), "visible"));
      await session.step(21, "And PhylocanvasGL viewer should fill its parent", () => fillsParent(page, el("PhylocanvasGL viewer")));
      await session.step(22, "And the table of PhylocanvasGL viewer should hold a tree with leaves \"leaf1, leaf2, leaf3\"", () => tableOfTreeLeaves(page, el("PhylocanvasGL viewer"), "leaf1, leaf2, leaf3"));
      await session.step(23, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on the file previews the same tree", async () => {
      await session.step(26, "When user clicks on \"Files > App Data > Dendrogram > data\" tree node inside browse tree", () => clickOn(page, el("\"Files > App Data > Dendrogram > data\" tree node inside browse tree")));
      await session.step(27, "And user clicks on nwk1.nwk link in gallery", () => clickOn(page, el("nwk1.nwk link in gallery")));
      await session.step(28, "Then PhylocanvasGL viewer should be visible", () => shouldBe(page, el("PhylocanvasGL viewer"), "visible"));
      await session.step(29, "And PhylocanvasGL viewer should fill its parent", () => fillsParent(page, el("PhylocanvasGL viewer")));
      await session.step(30, "And the table of PhylocanvasGL viewer should hold a tree with leaves \"leaf1, leaf2, leaf3\"", () => tableOfTreeLeaves(page, el("PhylocanvasGL viewer"), "leaf1, leaf2, leaf3"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The file handler makes a tree table with a Dendrogram of the same leaves", async () => {
      await session.step(34, "When user opens the newick file \"System:AppData/Dendrogram/data/nwk1.nwk\" with its file handler", () => openNewickFile(page, "System:AppData/Dendrogram/data/nwk1.nwk"));
      await session.step(35, "Then Dendrogram viewer should be visible", () => shouldBe(page, el("Dendrogram viewer"), "visible"));
      await session.step(36, "And table \"Table\" should have 5 rows", () => tableRows(page, "Table", 5));
      await session.step(37, "And table \"Table\" should have columns \"node, parent, leaf, distance\"", () => tableColumns(page, "Table", "node, parent, leaf, distance"));
      await session.step(38, "And the table should have tag \".newick\" equal to the text of \"System:AppData/Dendrogram/data/nwk1.nwk\" file", () => tableTagIsFile(page, ".newick", "System:AppData/Dendrogram/data/nwk1.nwk"));
      await session.step(39, "And the \"leaves\" reading of Dendrogram viewer should be \"leaf1, leaf2, leaf3\"", () => readingReads(page, "leaves", el("Dendrogram viewer"), "leaf1, leaf2, leaf3"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
