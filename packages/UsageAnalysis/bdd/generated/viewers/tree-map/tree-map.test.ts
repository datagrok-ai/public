/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree-map/tree-map.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tree-map]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaPainted, hasArea, hasNoArea, noErrors, painted, readingIs, readingReads, repainted, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree map splitting, nesting and the on-viewer selectors", () => {
  const session = feature(test, "features/viewers/tree-map/tree-map.feature", import.meta.url);
  test("Tree map splitting, nesting and the on-viewer selectors", {tag: ["@journey", "@viewers", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a tree map viewer", () => addViewer(page, "tree map"));
    await session.step(17, "Then the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(18, "And the \"error\" reading of tree map viewer should be \"\"", () => readingReads(page, "error", el("tree map viewer"), ""));
    await session.step(19, "And tree map viewer should be painted", () => painted(page, el("tree map viewer")));
    await run.scenario("The map opens split by the first column with 5 to 19 categories", async () => {
      await session.step(22, "Then the \"split columns\" reading of tree map viewer should be \"DIS_POP\"", () => readingReads(page, "split columns", el("tree map viewer"), "DIS_POP"));
      await session.step(23, "And the \"levels\" reading of tree map viewer should be 1", () => readingIs(page, "levels", el("tree map viewer"), 1));
      await session.step(24, "And the \"selectors\" reading of tree map viewer should be 2", () => readingIs(page, "selectors", el("tree map viewer"), 2));
      await session.step(25, "And the \"leaves\" reading of tree map viewer should be 6", () => readingIs(page, "leaves", el("tree map viewer"), 6));
      await session.step(26, "And the \"groups\" reading of tree map viewer should be 0", () => readingIs(page, "groups", el("tree map viewer"), 0));
      await session.step(27, "And the \"rows of RA\" reading of tree map viewer should be 434", () => readingIs(page, "rows of RA", el("tree map viewer"), 434));
      await session.step(28, "And the \"rows of Psoriasis\" reading of tree map viewer should be 204", () => readingIs(page, "rows of Psoriasis", el("tree map viewer"), 204));
      await session.step(29, "And the \"rows of PsA\" reading of tree map viewer should be 38", () => readingIs(page, "rows of PsA", el("tree map viewer"), 38));
      await session.step(30, "And tree map viewer should have a \"leaf RA\" area", () => hasArea(page, el("tree map viewer"), "leaf RA"));
      await session.step(31, "And tree map viewer should have a \"leaf body RA\" area", () => hasArea(page, el("tree map viewer"), "leaf body RA"));
      await session.step(32, "And the \"leaf RA\" area of tree map viewer should be painted", () => areaPainted(page, "leaf RA", el("tree map viewer")));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Splitting by RACE gives one rectangle per race holding exactly its rows", async () => {
      await session.step(36, "When user sets \"splitByColumnNames\" property of tree map viewer to \"RACE\"", () => setProperty(page, "splitByColumnNames", el("tree map viewer"), "RACE"));
      await session.step(37, "Then the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
      await session.step(38, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
      await session.step(39, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(40, "And the \"rows of Other\" reading of tree map viewer should be 62", () => readingIs(page, "rows of Other", el("tree map viewer"), 62));
      await session.step(41, "And the \"rows of Black\" reading of tree map viewer should be 27", () => readingIs(page, "rows of Black", el("tree map viewer"), 27));
      await session.step(42, "And the \"rows of Asian\" reading of tree map viewer should be 15", () => readingIs(page, "rows of Asian", el("tree map viewer"), 15));
      await session.step(43, "And the \"area of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "area of Caucasian", el("tree map viewer"), 896));
      await session.step(44, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(45, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(46, "And tree map viewer should have a \"leaf Caucasian\" area", () => hasArea(page, el("tree map viewer"), "leaf Caucasian"));
      await session.step(47, "And tree map viewer should not have a \"leaf RA\" area", () => hasNoArea(page, el("tree map viewer"), "leaf RA"));
      await session.step(48, "When user sets \"splitByColumnNames\" property of tree map viewer to \"DIS_POP\"", () => setProperty(page, "splitByColumnNames", el("tree map viewer"), "DIS_POP"));
      await session.step(49, "Then the \"leaves\" reading of tree map viewer should be 6", () => readingIs(page, "leaves", el("tree map viewer"), 6));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second level nests each race into its sexes and the group sums them", async () => {
      await session.step(53, "When user sets \"splitByColumnNames\" property of tree map viewer to \"RACE, SEX\"", () => setProperty(page, "splitByColumnNames", el("tree map viewer"), "RACE, SEX"));
      await session.step(54, "Then the \"levels\" reading of tree map viewer should be 2", () => readingIs(page, "levels", el("tree map viewer"), 2));
      await session.step(55, "And the \"selectors\" reading of tree map viewer should be 3", () => readingIs(page, "selectors", el("tree map viewer"), 3));
      await session.step(56, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
      await session.step(57, "And the \"groups\" reading of tree map viewer should be 4", () => readingIs(page, "groups", el("tree map viewer"), 4));
      await session.step(58, "And the \"rows of Caucasian | M\" reading of tree map viewer should be 416", () => readingIs(page, "rows of Caucasian | M", el("tree map viewer"), 416));
      await session.step(59, "And the \"rows of Caucasian | F\" reading of tree map viewer should be 480", () => readingIs(page, "rows of Caucasian | F", el("tree map viewer"), 480));
      await session.step(60, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(61, "And tree map viewer should have a \"group Caucasian\" area", () => hasArea(page, el("tree map viewer"), "group Caucasian"));
      await session.step(62, "And tree map viewer should have a \"leaf Caucasian | M\" area", () => hasArea(page, el("tree map viewer"), "leaf Caucasian | M"));
      await session.step(63, "And tree map viewer should not have a \"leaf Caucasian\" area", () => hasNoArea(page, el("tree map viewer"), "leaf Caucasian"));
      await session.step(64, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(65, "When user sets \"splitByColumnNames\" property of tree map viewer to \"DIS_POP\"", () => setProperty(page, "splitByColumnNames", el("tree map viewer"), "DIS_POP"));
      await session.step(66, "Then the \"levels\" reading of tree map viewer should be 1", () => readingIs(page, "levels", el("tree map viewer"), 1));
      await session.step(67, "And the \"groups\" reading of tree map viewer should be 0", () => readingIs(page, "groups", el("tree map viewer"), 0));
      await session.step(68, "And the \"selectors\" reading of tree map viewer should be 2", () => readingIs(page, "selectors", el("tree map viewer"), 2));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Column Selection Panel takes the selectors off the map", async () => {
      await session.step(72, "Then the \"selection panel shown\" reading of tree map viewer should be \"true\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "true"));
      await session.step(73, "And tree map viewer should have a \"split selector 1\" area", () => hasArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(74, "And tree map viewer should have a \"split selector 2\" area", () => hasArea(page, el("tree map viewer"), "split selector 2"));
      await session.step(75, "And tree map viewer should have a \"color selector\" area", () => hasArea(page, el("tree map viewer"), "color selector"));
      await session.step(76, "When user sets \"showColumnSelectionPanel\" property of tree map viewer to \"false\"", () => setProperty(page, "showColumnSelectionPanel", el("tree map viewer"), "false"));
      await session.step(77, "Then the \"selection panel shown\" reading of tree map viewer should be \"false\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "false"));
      await session.step(78, "And tree map viewer should not have a \"split selector 1\" area", () => hasNoArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(79, "And tree map viewer should not have a \"split selector 2\" area", () => hasNoArea(page, el("tree map viewer"), "split selector 2"));
      await session.step(80, "And tree map viewer should not have a \"color selector\" area", () => hasNoArea(page, el("tree map viewer"), "color selector"));
      await session.step(81, "And the \"levels\" reading of tree map viewer should be 1", () => readingIs(page, "levels", el("tree map viewer"), 1));
      await session.step(82, "When user sets \"showColumnSelectionPanel\" property of tree map viewer to \"true\"", () => setProperty(page, "showColumnSelectionPanel", el("tree map viewer"), "true"));
      await session.step(83, "Then the \"selection panel shown\" reading of tree map viewer should be \"true\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "true"));
      await session.step(84, "And tree map viewer should have a \"split selector 1\" area", () => hasArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the map", async () => {
      await session.step(88, "When user clicks on close icon of tree map viewer", () => clickOn(page, el("close icon of tree map viewer")));
      await session.step(89, "Then tree map viewer should be absent", () => shouldBe(page, el("tree map viewer"), "absent"));
      await session.step(90, "And the open tableview should have 0 tree map viewers", () => viewerCount(page, 0, "tree map"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
