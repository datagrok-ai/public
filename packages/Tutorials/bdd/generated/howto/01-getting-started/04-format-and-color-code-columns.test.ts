/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/04-format-and-color-code-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedAs, colorCodedCategorically} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/04-format-and-color-code-columns.feature", import.meta.url);
  test("Format and color-code columns", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "When user picks \"Format > scientific (1.60E2)\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Format > scientific (1.60E2)", "header HEIGHT", el("grid")), undefined, "Right-click a header, then Format: each choice shows a value in it");
    await session.step(21, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(22, "Then context menu should be absent", () => shouldBe(page, el("context menu"), "absent"));
    await session.step(23, "And \"HEIGHT\" column should have tag \"format\" equal to \"scientific\"", () => columnTag(page, "HEIGHT", "format", "scientific"));
    await session.step(25, "When user picks \"Format > max two digits after comma (73.2)\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Format > max two digits after comma (73.2)", "header WEIGHT", el("grid")), undefined, "At most two digits after the point");
    await session.step(26, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(27, "Then context menu should be absent", () => shouldBe(page, el("context menu"), "absent"));
    await session.step(28, "And \"WEIGHT\" column should have tag \"format\" equal to \"max two digits after comma\"", () => columnTag(page, "WEIGHT", "format", "max two digits after comma"));
    await session.step(30, "When user picks \"Format > dd.MM.yyyy (11.09.1990)\" from the context menu of the \"header STARTED\" area of grid", () => pickFromAreaContextMenu(page, "Format > dd.MM.yyyy (11.09.1990)", "header STARTED", el("grid")), undefined, "A date column offers date formats");
    await session.step(31, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(32, "Then context menu should be absent", () => shouldBe(page, el("context menu"), "absent"));
    await session.step(33, "And \"STARTED\" column should have tag \"format\" equal to \"dd.MM.yyyy\"", () => columnTag(page, "STARTED", "format", "dd.MM.yyyy"));
    await session.step(35, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")), undefined, "Color Coding > Linear: shades from low to high");
    await session.step(36, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(37, "Then context menu should be absent", () => shouldBe(page, el("context menu"), "absent"));
    await session.step(38, "And \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
    await session.step(40, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header RACE", el("grid")), undefined, "Color Coding > Categorical: a color for each category");
    await session.step(41, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(42, "Then context menu should be absent", () => shouldBe(page, el("context menu"), "absent"));
    await session.step(43, "And \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
});
