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
import {addWideViewer, hideGridColumn, hoverMore, hoverShortened, legendCorner, legendDrawn, legendNotDrawn, pointerOffLegend, tooltipListsCurves, tooltipWholeName} from '../../bindings/curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, expand, selectIn, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickOtherCell, pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The legend of a dose-response chart names the curves without hiding them", () => {
  const session = feature(test, "features/guides/curves-legend.feature", import.meta.url);
  test("Name the curves, read a shortened name and the folded rows, then hide and show the legend", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(20, "And the grid of the current view does not show the \"Molecule\" column", () => hideGridColumn(page, "Molecule"));
    await session.step(21, "And user selects rows where \"Assay Name\" is \"Collapse Assay\"", () => selectWhereIs(page, "Assay Name", "Collapse Assay"));
    await session.step(22, "And user adds a MultiCurveViewer viewer of the selected rows taking 45 percent of the view width", () => addWideViewer(page, 45));
    await session.step(23, "When user clicks on a \"AUC\" cell of grid other than the current one", () => clickOtherCell(page, "AUC", el("grid")));
    await session.step(24, "And user picks \"Properties...\" from the viewer menu of MultiCurveViewer viewer", () => pickFromViewerMenu(page, "Properties...", el("MultiCurveViewer viewer")));
    await session.step(25, "And user selects \"Compound|Assay|Target\" in \"Legend\" property in context panel", () => selectIn(page, "Compound|Assay|Target", el("\"Legend\" property in context panel")));
    await session.step(26, "Then \"Legend\" property of MultiCurveViewer viewer should be \"Compound|Assay|Target\"", () => propertyShouldBe(page, "Legend", el("MultiCurveViewer viewer"), "Compound|Assay|Target"));
    await session.step(27, "And the curves legend of MultiCurveViewer viewer should be in the bottom left corner", () => legendCorner(page, el("MultiCurveViewer viewer"), "bottom", "left"));
    await session.step(28, "When user hovers over a shortened name in the curves legend of MultiCurveViewer viewer", () => hoverShortened(page, el("MultiCurveViewer viewer")));
    await session.step(29, "Then the tooltip should show the whole of the shortened name", () => tooltipWholeName(page));
    await session.step(30, "When user hovers over the \"+N more\" row of the curves legend of MultiCurveViewer viewer", () => hoverMore(page, el("MultiCurveViewer viewer")));
    await session.step(31, "Then the tooltip should list the curves of the legend, each in its own color", () => tooltipListsCurves(page));
    await session.step(32, "When user moves the pointer off the curves legend of MultiCurveViewer viewer", () => pointerOffLegend(page, el("MultiCurveViewer viewer")));
    await session.step(33, "And user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")));
    await session.step(34, "And user unchecks \"Show Legend\" property in context panel", () => uncheck(page, el("\"Show Legend\" property in context panel")));
    await session.step(35, "Then the curves legend of MultiCurveViewer viewer should not be drawn", () => legendNotDrawn(page, el("MultiCurveViewer viewer")));
    await session.step(36, "When user checks \"Show Legend\" property in context panel", () => check(page, el("\"Show Legend\" property in context panel")));
    await session.step(37, "Then the curves legend of MultiCurveViewer viewer should be drawn", () => legendDrawn(page, el("MultiCurveViewer viewer")));
  });
});
