/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/scatter-plot-marker-labels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, readingAtLeast} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Label scatter plot markers with the structure, the ID and a potency value", () => {
  const session = feature(test, "features/guides/scatter-plot-marker-labels.feature", import.meta.url);
  test("Pick the label columns of a scatter plot", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(13, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","Chemical Space X"],["Y","Chemical Space Y"]]), [["X","Chemical Space X"],["Y","Chemical Space Y"]]);
    await session.step(16, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(17, "Given \"Labels\" category in context panel is expanded", () => isExpanded(page, el("\"Labels\" category in context panel")));
    await session.step(18, "When user clicks on \"...\" button in \"Label Columns\" property in context panel", () => clickOn(page, el("\"...\" button in \"Label Columns\" property in context panel")));
    await session.step(19, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(20, "When user toggles the \"Structure\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Structure", el("\"Select columns...\" dialog")));
    await session.step(21, "And user toggles the \"Id\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Id", el("\"Select columns...\" dialog")));
    await session.step(22, "And user types \"Cellular assay 1\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "Cellular assay 1", el("\"Search\" input in \"Select columns...\" dialog")));
    await session.step(23, "And user toggles the \"Cellular assay 1\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Cellular assay 1", el("\"Select columns...\" dialog")));
    await session.step(24, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(25, "Then the \"labels shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 1));
  });
});
