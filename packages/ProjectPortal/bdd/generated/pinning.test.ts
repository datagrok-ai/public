/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/pinning.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addressEndsWith, notPinned, openApp, openHome, openTile} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Add a program to favorites and find it in the Workspace and on Home", () => {
  const session = feature(test, "features/pinning.feature", import.meta.url);
  test("Add a discovery project to favorites, see it on Home, open it from there, remove it", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given the user has not pinned \"EXS-KIN1\" in the \"EX-SMD\" portfolio", () => notPinned(page, "EXS-KIN1", "EX-SMD"));
    await session.step(12, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(13, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(14, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(15, "Then pin button should contain text \"Add to favorites\"", () => shouldContainText(page, el("pin button"), "Add to favorites"));
    await session.step(16, "When user clicks on pin button", () => clickOn(page, el("pin button")));
    await session.step(17, "Then pin button should contain text \"Remove from favorites\"", () => shouldContainText(page, el("pin button"), "Remove from favorites"));
    await session.step(18, "When user clicks on board link", () => clickOn(page, el("board link")));
    await session.step(19, "Then EXS-KIN1 pinned mark should be visible", () => shouldBe(page, el("EXS-KIN1 pinned mark"), "visible"));
    await session.step(20, "When user opens Home", () => openHome(page));
    await session.step(21, "Then my programs widget should contain text \"EXS-KIN1\"", () => shouldContainText(page, el("my programs widget"), "EXS-KIN1"));
    await session.step(22, "And EXS-KIN1 tile should contain text \"Example Small Molecule Discovery\"", () => shouldContainText(page, el("EXS-KIN1 tile"), "Example Small Molecule Discovery"));
    await session.step(23, "When user clicks on EXS-KIN1 tile", () => clickOn(page, el("EXS-KIN1 tile")));
    await session.step(24, "Then program panel should contain text \"Example Small Molecule Discovery\"", () => shouldContainText(page, el("program panel"), "Example Small Molecule Discovery"));
    await session.step(25, "And pin button should contain text \"Remove from favorites\"", () => shouldContainText(page, el("pin button"), "Remove from favorites"));
    await session.step(26, "When user double-clicks on the \"EXS-KIN1\" tile", () => openTile(page, "EXS-KIN1"));
    await session.step(27, "Then the address should end with \"/Portal/EX-SMD/EXS-KIN1\"", () => addressEndsWith(page, "/Portal/EX-SMD/EXS-KIN1"));
    await session.step(28, "And pin button should contain text \"Remove from favorites\"", () => shouldContainText(page, el("pin button"), "Remove from favorites"));
    await session.step(29, "When user clicks on pin button", () => clickOn(page, el("pin button")));
    await session.step(30, "Then pin button should contain text \"Add to favorites\"", () => shouldContainText(page, el("pin button"), "Add to favorites"));
    await session.step(31, "When user clicks on board link", () => clickOn(page, el("board link")));
    await session.step(32, "Then EXS-KIN1 pinned mark should not be visible", () => shouldNotBe(page, el("EXS-KIN1 pinned mark"), "visible"));
  });
});
