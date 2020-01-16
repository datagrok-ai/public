const utils = require("./utils.js");
const tests = {};
tests["console-test"] = async (driver, vars, opts = {}) => {
  await driver.get((new URL(`/datasets?selenium=true&q=`, BASE_URL)).href);
  await driver.wait(until.elementLocated(By.xpath(`//div[@id=\'signup-login-fields\']/div/div/input`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@id=\'signup-login-fields\']/div/div/input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`selenium`);
    });
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@id=\'signup-login-fields\']/div/div[2]/input`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@id=\'signup-login-fields\']/div/div[2]/input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`selenium`);
    });
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@id=\'signup-login-fields\']/div[2]/button`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@id=\'signup-login-fields\']/div[2]/button`)).then(element => {
    return element.click();
  });
  await driver.sleep(4000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Tools\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Tools\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Tools---Console\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Tools---Console\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//pre[contains(.,\'CmdContextHelp()\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'CmdContextHelp()\')]`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//pre[contains(.,\'CmdConsole()\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'CmdConsole()\')]`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-trash-alt\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-trash-alt\']`)).then(element => {
    return element.click();
  });
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'CmdContextHelp()\')]`))).resolves.not.toBePresent();
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'CmdConsole()\')]`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.name(`div-SeleniumTestDataP1`)), configuration.timeout);
  await driver.findElement(By.name(`div-SeleniumTestDataP1`)).then(element => {
    return driver.actions({
      bridge: true
    }).doubleClick(element).perform();
  });
  await driver.sleep(2000);
  await driver.wait(until.elementLocated(By.xpath(`//pre[contains(.,\'  result: demog (5850 rows, 10 columns)\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'  result: demog (5850 rows, 10 columns)\')]`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//pre[contains(.,\'  result: Project \"Selenium test data p1\"\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'  result: Project \"Selenium test data p1\"\')]`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`1+1.5`);
    });
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.wait(until.elementLocated(By.xpath(`//pre[contains(.,\'  result: 2.5\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'  result: 2.5\')]`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`a=10`);
    });
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.sleep(1600);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`b=\"selenium-test\"`);
    });
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`a`);
    });
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.wait(until.elementLocated(By.xpath(`//pre[contains(.,\'  value: 10\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'  value: 10\')]`))).resolves.toBePresent();
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`b`);
    });
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.wait(until.elementLocated(By.xpath(`//pre[contains(.,\'  value: selenium-test\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//pre[contains(.,\'  value: selenium-test\')]`))).resolves.toBePresent();
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`SelectRows(, HasNulls(Named([\"USUBJID\", \"AGE\", \"SEX\", \"RACE\", \"DIS_POP\", \"HEIGHT\", \"WEIGHT\", \"DEMOG\", \"CONTROL\", \"STARTED\"])))`);
    });
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.name(`span-selected`)), configuration.timeout);
  await expect(driver.findElement(By.name(`span-selected`))).resolves.toHaveText(`753 selected rows`);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-list\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-list\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`.panel-base:nth-child(3) > .panel-titlebar`)), configuration.timeout);
  await expect(driver.findElement(By.css(`.panel-base:nth-child(3) > .panel-titlebar`))).resolves.toHaveText(`Variables`);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.css(`canvas:nth-child(4)`)), configuration.timeout);
  await driver.findElement(By.css(`canvas:nth-child(4)`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`pc`);
    });
  });
  await driver.sleep(1000);
  await driver.executeScript(`let element = document.querySelector('.d4-console-input > input'); let event = new KeyboardEvent('keydown', {keyCode: 9}); element.dispatchEvent(event);`);
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.sleep(2000);
  await driver.wait(until.elementLocated(By.name(`viewer-Pie-chart`)), configuration.timeout);
  await expect(driver.findElements(By.name(`viewer-Pie-chart`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`CmdC`);
    });
  });
  await driver.executeScript(`let element = document.querySelector('.d4-console-input > input'); let event = new KeyboardEvent('keydown', {keyCode: 9}); element.dispatchEvent(event);`);
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.css(`.d4-console-input > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-console-input > input`)).then(element => {
    return element.sendKeys(Key["ENTER"]);
  });
  await driver.sleep(2000);
  await driver.wait(until.elementLocated(By.css(`.d4-dialog`)), configuration.timeout);
  await expect(driver.findElements(By.css(`.d4-dialog`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.css(`.d4-dialog-header`)), configuration.timeout);
  await expect(driver.findElement(By.css(`.d4-dialog-header`))).resolves.toHaveText(`Cluster Data`);
  await driver.sleep(1000);
  await driver.executeScript(`let element = window; let event = new KeyboardEvent('keydown', {keyCode: 192}); element.dispatchEvent(event);`);
  await driver.sleep(1000);
  await expect(driver.findElements(By.css(`.d4-console-wrapper`))).resolves.not.toBePresent();
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
}
module.exports = tests;