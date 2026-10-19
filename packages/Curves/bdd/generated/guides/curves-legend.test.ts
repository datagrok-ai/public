/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/curves-legend.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addWideViewer} from '../../bindings/curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, expand, selectIn, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {removeColumn, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, hoverFirstArea, pointerAway, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickOtherCell, pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The legend of a dose-response chart names the curves without hiding them", () => {
  const session = feature(test, "features/guides/curves-legend.feature", import.meta.url);
  test("Name the curves, read a shortened name and the folded rows, then hide and show the legend", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(25, "And user removes \"Molecule\" column", () => removeColumn(page, "Molecule"));
    await session.step(26, "And user selects rows where \"Assay Name\" is \"Collapse Assay\"", () => selectWhereIs(page, "Assay Name", "Collapse Assay"));
    await session.step(27, "And user adds a MultiCurveViewer viewer of the selected rows taking 45 percent of the view width", () => addWideViewer(page, 45));
    await session.step(28, "When user clicks on a \"AUC\" cell of grid other than the current one", () => clickOtherCell(page, "AUC", el("grid")));
    await session.step(29, "And user picks \"Properties...\" from the viewer menu of MultiCurveViewer viewer", () => pickFromViewerMenu(page, "Properties...", el("MultiCurveViewer viewer")));
    await session.step(30, "And user selects \"Compound|Assay|Target\" in \"Legend\" property in context panel", () => selectIn(page, "Compound|Assay|Target", el("\"Legend\" property in context panel")));
    await session.step(31, "Then the \"legend corner\" reading of MultiCurveViewer viewer should be \"bottom left\"", () => readingReads(page, "legend corner", el("MultiCurveViewer viewer"), "bottom left"));
    await session.step(32, "When user hovers over the first \"shortened legend row\" area of MultiCurveViewer viewer", () => hoverFirstArea(page, "shortened legend row", el("MultiCurveViewer viewer")));
    await session.step(33, "Then tooltip should contain text \"| Collapse Assay | target_10\"", () => shouldContainText(page, el("tooltip"), "| Collapse Assay | target_10"));
    await session.step(34, "When user hovers over the \"legend more\" area of MultiCurveViewer viewer", () => hoverArea(page, "legend more", el("MultiCurveViewer viewer")));
    await session.step(35, "Then tooltip should contain text \"compound_id 20 | Collapse Assay | target_10\"", () => shouldContainText(page, el("tooltip"), "compound_id 20 | Collapse Assay | target_10"));
    await session.step(36, "When user moves the pointer away from MultiCurveViewer viewer", () => pointerAway(page, el("MultiCurveViewer viewer")));
    await session.step(37, "And user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")));
    await session.step(38, "Then \"Show Legend\" property in context panel should be checked", () => shouldBe(page, el("\"Show Legend\" property in context panel"), "checked"));
    await session.step(39, "When user unchecks \"Show Legend\" property in context panel", () => uncheck(page, el("\"Show Legend\" property in context panel")));
    await session.step(40, "Then the \"legend shown\" reading of MultiCurveViewer viewer should be \"false\"", () => readingReads(page, "legend shown", el("MultiCurveViewer viewer"), "false"));
    await session.step(41, "When user checks \"Show Legend\" property in context panel", () => check(page, el("\"Show Legend\" property in context panel")));
    await session.step(42, "Then the \"legend shown\" reading of MultiCurveViewer viewer should be \"true\"", () => readingReads(page, "legend shown", el("MultiCurveViewer viewer"), "true"));
  });
});
