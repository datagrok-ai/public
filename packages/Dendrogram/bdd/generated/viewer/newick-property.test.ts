/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewer/newick-property.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [dendrogram.viewer.newick-property]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorLinear, setTableTag} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, moreInk, noErrors, painted, paintedInColors, readingReads, repainted, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Dendrogram viewer drawn from a newick it is given", () => {
  const session = feature(test, "features/viewer/newick-property.feature", import.meta.url);
  test("The Dendrogram viewer drawn from a newick it is given", {tag: ["@journey", "@realizes:dendrogram.viewer.newick-property", "@known-failure:52"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And user opens a table \"leaves4\" with:", () => openTableOf(page, "leaves4", [["leaf","value"],["a","1"],["b","2"],["c","3"],["d","4"]]), [["leaf","value"],["a","1"],["b","2"],["c","3"],["d","4"]]);
    await session.step(35, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(36, "And user sets the \".newick\" tag of the table to \"((d:1,c:1):1,(b:1,a:1):1);\"", () => setTableTag(page, ".newick", "((d:1,c:1):1,(b:1,a:1):1);"));
    await session.step(37, "And user sets the \".newick-alt\" tag of the table to \"(((a:1,b:1):1,c:1):1,d:1);\"", () => setTableTag(page, ".newick-alt", "(((a:1,b:1):1,c:1):1,d:1);"));
    await session.step(38, "And user adds a Dendrogram viewer with:", () => addViewerWith(page, "Dendrogram", [["newick","((a:1,b:1):1,(c:1,d:1):1);"],["nodeColumnName","leaf"]]), [["newick","((a:1,b:1):1,(c:1,d:1):1);"],["nodeColumnName","leaf"]]);
    await run.scenario("The newick the viewer is given wins over the table's tags", async () => {
      await session.step(43, "Then Dendrogram viewer should be visible", () => shouldBe(page, el("Dendrogram viewer"), "visible"));
      await session.step(44, "And the \"newick\" reading of Dendrogram viewer should be \"((a:1,b:1):1,(c:1,d:1):1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "((a:1,b:1):1,(c:1,d:1):1);"));
      await session.step(45, "And the \"leaves\" reading of Dendrogram viewer should be \"a, b, c, d\"", () => readingReads(page, "leaves", el("Dendrogram viewer"), "a, b, c, d"));
      await session.step(46, "And Dendrogram viewer should be painted", () => painted(page, el("Dendrogram viewer")));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Labels draws the leaf names", async () => {
      await session.step(51, "When user sets \"showLabels\" property of Dendrogram viewer to \"true\"", () => setProperty(page, "showLabels", el("Dendrogram viewer"), "true"));
      await session.step(52, "Then Dendrogram viewer should have more ink than before", () => moreInk(page, el("Dendrogram viewer")));
      await session.step(53, "When user sets \"showLabels\" property of Dendrogram viewer to \"false\"", () => setProperty(page, "showLabels", el("Dendrogram viewer"), "false"));
    }, {knownFailure: true, line: 52});
    await run.scenario("Without its own newick the viewer draws the tag Newick Tag names, else .newick", async () => {
      await session.step(56, "When user sets \"newick\" property of Dendrogram viewer to \"\"", () => setProperty(page, "newick", el("Dendrogram viewer"), ""));
      await session.step(57, "Then the \"newick\" reading of Dendrogram viewer should be \"((d:1,c:1):1,(b:1,a:1):1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "((d:1,c:1):1,(b:1,a:1):1);"));
      await session.step(58, "And the \"leaves\" reading of Dendrogram viewer should be \"d, c, b, a\"", () => readingReads(page, "leaves", el("Dendrogram viewer"), "d, c, b, a"));
      await session.step(59, "When user clicks on settings icon of Dendrogram viewer", () => clickOn(page, el("settings icon of Dendrogram viewer")));
      await session.step(60, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(61, "When user selects \".newick-alt\" in \"Newick Tag\" property in context panel", () => selectIn(page, ".newick-alt", el("\"Newick Tag\" property in context panel")));
      await session.step(62, "Then the \"newick\" reading of Dendrogram viewer should be \"(((a:1,b:1):1,c:1):1,d:1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "(((a:1,b:1):1,c:1):1,d:1);"));
      await session.step(63, "And the \"leaves\" reading of Dendrogram viewer should be \"a, b, c, d\"", () => readingReads(page, "leaves", el("Dendrogram viewer"), "a, b, c, d"));
      await session.step(64, "When user sets \"newickTag\" property of Dendrogram viewer to \"\"", () => setProperty(page, "newickTag", el("Dendrogram viewer"), ""));
      await session.step(65, "And user sets \"newick\" property of Dendrogram viewer to \"((a:1,b:1):1,(c:1,d:1):1);\"", () => setProperty(page, "newick", el("Dendrogram viewer"), "((a:1,b:1):1,(c:1,d:1):1);"));
      await session.step(66, "Then the \"newick\" reading of Dendrogram viewer should be \"((a:1,b:1):1,(c:1,d:1):1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "((a:1,b:1):1,(c:1,d:1):1);"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A color column colors the leaves by its color coding", async () => {
      await session.step(70, "When user colors \"value\" column linearly from \"#0000ff\" to \"#ff0000\"", () => colorLinear(page, "value", "#0000ff", "#ff0000"));
      await session.step(71, "And user sets \"colorColumnName\" property of Dendrogram viewer to \"value\"", () => setProperty(page, "colorColumnName", el("Dendrogram viewer"), "value"));
      await session.step(72, "Then Dendrogram viewer should have repainted", () => repainted(page, el("Dendrogram viewer")));
      await session.step(73, "And Dendrogram viewer should be painted in at least 3 colors", () => paintedInColors(page, el("Dendrogram viewer"), 3));
      await session.step(74, "When user sets \"colorColumnName\" property of Dendrogram viewer to \"\"", () => setProperty(page, "colorColumnName", el("Dendrogram viewer"), ""));
      await session.step(75, "Then Dendrogram viewer should have repainted", () => repainted(page, el("Dendrogram viewer")));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The style properties restyle the tree", async () => {
      await session.step(79, "When user sets \"lineWidth\" property of Dendrogram viewer to \"8\"", () => setProperty(page, "lineWidth", el("Dendrogram viewer"), "8"));
      await session.step(80, "Then Dendrogram viewer should have repainted", () => repainted(page, el("Dendrogram viewer")));
      await session.step(81, "When user sets \"nodeSize\" property of Dendrogram viewer to \"12\"", () => setProperty(page, "nodeSize", el("Dendrogram viewer"), "12"));
      await session.step(82, "Then Dendrogram viewer should have repainted", () => repainted(page, el("Dendrogram viewer")));
      await session.step(83, "When user sets \"showGrid\" property of Dendrogram viewer to \"true\"", () => setProperty(page, "showGrid", el("Dendrogram viewer"), "true"));
      await session.step(84, "Then Dendrogram viewer should have repainted", () => repainted(page, el("Dendrogram viewer")));
      await session.step(85, "When user sets \"mainColor\" property of Dendrogram viewer to \"#ff0000\"", () => setProperty(page, "mainColor", el("Dendrogram viewer"), "#ff0000"));
      await session.step(86, "Then Dendrogram viewer should have repainted", () => repainted(page, el("Dendrogram viewer")));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
