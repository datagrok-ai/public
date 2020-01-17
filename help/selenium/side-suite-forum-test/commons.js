const utils = require("./utils.js");
const tests = {};
tests["forum-test"] = async (driver, vars, opts = {}) => {
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
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Help\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Help\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Help---Forum\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Help---Forum\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.id(`community-forum`)), configuration.timeout);
  await expect(driver.findElements(By.id(`community-forum`))).resolves.toBePresent();
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.name(`button-NEW-TOPIC`)), configuration.timeout);
  await driver.findElement(By.name(`button-NEW-TOPIC`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.id(`subject`)), configuration.timeout);
  await driver.findElement(By.id(`subject`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.id(`subject`)), configuration.timeout);
  await driver.findElement(By.id(`subject`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`Topic by Selenium`);
    });
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.css(`.grok-comments-post-input`)), configuration.timeout);
  await driver.findElement(By.css(`.grok-comments-post-input`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`.grok-comments-post-input`)), configuration.timeout);
  await driver.findElement(By.css(`.grok-comments-post-input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`Text by Selenium`);
    });
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.name(`button-CREATE-TOPIC`)), configuration.timeout);
  await driver.findElement(By.name(`button-CREATE-TOPIC`)).then(element => {
    return element.click();
  });
  await driver.sleep(2000);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//td[contains(.,\'Topic by Selenium\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//td[contains(.,\'Topic by Selenium\')]`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-icon-user\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-icon-user\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.xpath(`//label[@name=\'label-Chats\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//label[@name=\'label-Chats\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//label[@name=\'label-Topic-by-Selenium\']`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//label[@name=\'label-Topic-by-Selenium\']`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Unwatch\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Unwatch\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Delete\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Delete\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await expect(driver.findElements(By.xpath(`//label[@name=\'label-Topic by Selenium\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Help\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Help\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Help---Forum\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Help---Forum\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await expect(driver.findElements(By.xpath(`//td[contains(.,\'Topic by Selenium\')]`))).resolves.not.toBePresent();
}
module.exports = tests;