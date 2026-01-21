var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import { testData } from './dataframe-utils';
import { changeOptionsSaveLayout, filterAsync, loadLayout, selectFilterChangeCurrent, testViewerInternal } from './test-viewer-utils';
const STANDART_TIMEOUT = 30000;
const BENCHMARK_TIMEOUT = 10800000;
const stdLog = console.log.bind(console);
const stdInfo = console.info.bind(console);
const stdWarn = console.warn.bind(console);
const stdError = console.error.bind(console);
export const tests = {};
const autoTestsCatName = 'Auto Tests';
const demoCatName = 'Demo';
const detectorsCatName = 'Detectors';
const coreCatName = 'Core';
const wasRegistered = {};
export let currentCategory;
export var assure;
(function (assure) {
    function notNull(value, name) {
        if (value == null)
            throw new Error(`${name == null ? 'Value' : name} not defined`);
    }
    assure.notNull = notNull;
})(assure || (assure = {}));
export class TestContext {
    constructor(catchUnhandled, report, returnOnFail) {
        this.catchUnhandled = true;
        this.report = false;
        this.returnOnFail = false;
        if (catchUnhandled !== undefined)
            this.catchUnhandled = catchUnhandled;
        if (report !== undefined)
            this.report = report;
        if (returnOnFail !== undefined)
            this.returnOnFail = returnOnFail;
    }
    ;
}
export class Test {
    constructor(category, name, test, options) {
        var _a;
        this.category = category;
        this.name = name;
        options !== null && options !== void 0 ? options : (options = {});
        (_a = options.timeout) !== null && _a !== void 0 ? _a : (options.timeout = STANDART_TIMEOUT);
        this.options = options;
        this.test = () => __awaiter(this, void 0, void 0, function* () {
            return new Promise((resolve, reject) => __awaiter(this, void 0, void 0, function* () {
                var _b;
                let result = '';
                try {
                    if (DG.Test.isInDebug)
                        debugger;
                    let res = yield test();
                    try {
                        result = (_b = res === null || res === void 0 ? void 0 : res.toString()) !== null && _b !== void 0 ? _b : '';
                    }
                    catch (e) {
                        result = 'Can\'t convert test\'s result to string';
                        console.error(`Can\'t convert test\'s result to string in the ${this.category}:${this.name} test`);
                    }
                }
                catch (e) {
                    reject(e);
                }
                resolve(result);
            }));
        });
    }
}
export class Category {
}
export class NodeTestExecutionOptions {
}
export class TestExecutionOptions {
}
export function testEvent(event, handler, trigger, ms = 0, reason = `timeout`) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => {
            const sub = event.subscribe((args) => {
                try {
                    handler(args);
                    resolve('OK');
                }
                catch (e) {
                    reject(e);
                }
                finally {
                    sub.unsubscribe();
                    clearTimeout(timeout);
                }
            });
            const timeout = setTimeout(() => {
                sub.unsubscribe();
                // eslint-disable-next-line prefer-promise-reject-errors
                reject(reason);
            }, ms);
            trigger();
        });
    });
}
export function testEventAsync(event, handler, trigger, ms = 0, reason = `timeout`) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => {
            const sub = event.subscribe((args) => {
                handler(args).then(() => {
                    resolve('OK');
                }).catch((e) => {
                    reject(e);
                }).finally(() => {
                    sub.unsubscribe();
                    clearTimeout(timeout);
                });
            });
            const timeout = setTimeout(() => {
                sub.unsubscribe();
                // eslint-disable-next-line prefer-promise-reject-errors
                reject(reason);
            }, ms);
            trigger();
        });
    });
}
export function test(name, test, options) {
    if (tests[currentCategory] == undefined)
        tests[currentCategory] = {};
    if (tests[currentCategory].tests == undefined)
        tests[currentCategory].tests = [];
    tests[currentCategory].tests.push(new Test(currentCategory, name, test, options));
}
/* Tests two objects for equality, throws an exception if they are not equal. */
export function expect(actual, expected = true, error) {
    if (error)
        error = `${error}, `;
    else
        error = '';
    if (actual !== expected)
        throw new Error(`${error}Expected "${expected}", got "${actual}"`);
}
export function expectFloat(actual, expected, tolerance = 0.001, error) {
    if ((actual === Number.POSITIVE_INFINITY && expected === Number.POSITIVE_INFINITY) ||
        (actual === Number.NEGATIVE_INFINITY && expected === Number.NEGATIVE_INFINITY) ||
        (actual === Number.NaN && expected === Number.NaN) || (isNaN(actual) && isNaN(expected)))
        return;
    const areEqual = Math.abs(actual - expected) < tolerance;
    expect(areEqual, true, `${error !== null && error !== void 0 ? error : ''} (tolerance = ${tolerance})`);
    if (!areEqual)
        throw new Error(`Expected ${expected}, got ${actual} (tolerance = ${tolerance})`);
}
export function expectTable(actual, expected, error) {
    const expectedRowCount = expected.rowCount;
    const actualRowCount = actual.rowCount;
    expect(actualRowCount, expectedRowCount, `${error !== null && error !== void 0 ? error : ''}, row count`);
    for (const column of expected.columns) {
        const actualColumn = actual.columns.byName(column.name);
        if (actualColumn == null)
            throw new Error(`Column ${column.name} not found`);
        if (actualColumn.type != column.type)
            throw new Error(`Column ${column.name} type expected ${column.type} got ${actualColumn.type}`);
        for (let i = 0; i < expectedRowCount; i++) {
            const value = column.get(i);
            const actualValue = actualColumn.get(i);
            if (column.type == DG.TYPE.FLOAT)
                expectFloat(actualValue, value, 0.0001, error);
            else if (column.type == DG.TYPE.DATE_TIME)
                expect(actualValue.isSame(value), true, error);
            else
                expect(actualValue, value, error);
        }
    }
}
export function expectObject(actual, expected) {
    for (const [expectedKey, expectedValue] of Object.entries(expected)) {
        if (!actual.hasOwnProperty(expectedKey))
            throw new Error(`Expected property "${expectedKey}" not found`);
        const actualValue = actual[expectedKey];
        if (actualValue instanceof Array && expectedValue instanceof Array)
            expectArray(actualValue, expectedValue);
        else if (actualValue instanceof Object && expectedValue instanceof Object)
            expectObject(actualValue, expectedValue);
        else if (Number.isFinite(actualValue) && Number.isFinite(expectedValue))
            expectFloat(actualValue, expectedValue);
        else if (actualValue != expectedValue)
            throw new Error(`Expected (${expectedValue}) for key '${expectedKey}', got (${actualValue})`);
    }
}
export function expectArray(actual, expected) {
    const actualLength = actual.length;
    const expectedLength = expected.length;
    if (actualLength != expectedLength) {
        throw new Error(`Arrays are of different length: actual array length is ${actualLength} ` +
            `and expected array length is ${expectedLength}`);
    }
    for (let i = 0; i < actualLength; i++) {
        if (actual[i] instanceof Array && expected[i] instanceof Array)
            expectArray(actual[i], expected[i]);
        else if (actual[i] instanceof Object && expected[i] instanceof Object)
            expectObject(actual[i], expected[i]);
        else if (actual[i] != expected[i])
            throw new Error(`Expected ${expected[i]} at position ${i}, got ${actual[i]}`);
    }
}
/* Defines a test suite. */
export function category(category, tests_, options) {
    var _a;
    currentCategory = category;
    tests_();
    if (tests[currentCategory]) {
        tests[currentCategory].clear = (_a = options === null || options === void 0 ? void 0 : options.clear) !== null && _a !== void 0 ? _a : true;
        tests[currentCategory].timeout = options === null || options === void 0 ? void 0 : options.timeout;
        tests[currentCategory].benchmarks = options === null || options === void 0 ? void 0 : options.benchmarks;
        tests[currentCategory].stressTests = options === null || options === void 0 ? void 0 : options.stressTests;
        tests[currentCategory].owner = options === null || options === void 0 ? void 0 : options.owner;
    }
}
/* Defines a function to be executed before the tests in this category are executed. */
export function before(before) {
    if (tests[currentCategory] == undefined)
        tests[currentCategory] = {};
    tests[currentCategory].before = before;
}
/* Defines a function to be executed after the tests in this category are executed. */
export function after(after) {
    if (tests[currentCategory] == undefined)
        tests[currentCategory] = {};
    tests[currentCategory].after = after;
}
function addNamespace(s, f) {
    return s.replace(new RegExp(f.name, 'gi'), f.nqName);
}
export function initAutoTests(package_, module) {
    var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k;
    return __awaiter(this, void 0, void 0, function* () {
        const packageId = package_.id;
        if (wasRegistered[packageId])
            return;
        const moduleTests = module ? module.tests : tests;
        if (package_.name === 'DevTools' || (!!module && module._package.name === 'DevTools')) {
            for (const f of window.dartTests) {
                const arr = f.name.split(/\s*\|\s*!/g);
                let name = (_a = arr.pop()) !== null && _a !== void 0 ? _a : f.name;
                let cat = arr.length ? coreCatName + ': ' + arr.join(': ') : coreCatName;
                let fullName = name.split(' | ');
                name = fullName[fullName.length - 1];
                fullName.unshift(cat);
                fullName.pop();
                cat = fullName.join(': ');
                if (moduleTests[cat] === undefined)
                    moduleTests[cat] = { tests: [], clear: true };
                moduleTests[cat].tests.push(new Test(cat, name, f.test, { isAggregated: false, timeout: (_c = (_b = f.options) === null || _b === void 0 ? void 0 : _b.timeout) !== null && _c !== void 0 ? _c : STANDART_TIMEOUT, skipReason: (_d = f.options) === null || _d === void 0 ? void 0 : _d.skipReason, owner: (_e = f.options) === null || _e === void 0 ? void 0 : _e.owner, benchmark: (_g = (_f = f.options) === null || _f === void 0 ? void 0 : _f.benchmark) !== null && _g !== void 0 ? _g : false }));
            }
        }
        const moduleAutoTests = [];
        const moduleDemo = [];
        const moduleDetectors = [];
        const packFunctions = yield grok.dapi.functions.filter(`package.id = "${packageId}"`).list();
        const reg = new RegExp(/skip:\s*([^,\s]+)|wait:\s*(\d+)|cat:\s*([^,\s]+)|timeout:\s*(\d+)/g);
        for (const f of packFunctions) {
            const tests = f.options['test'];
            const demo = f.options['demoPath'];
            if ((tests && Array.isArray(tests) && tests.length)) {
                for (let i = 0; i < tests.length; i++) {
                    const res = tests[i].matchAll(reg);
                    const map = {};
                    Array.from(res).forEach((arr) => {
                        if (arr[0].startsWith('skip'))
                            map['skip'] = arr[1];
                        else if (arr[0].startsWith('wait'))
                            map['wait'] = parseInt(arr[2]);
                        else if (arr[0].startsWith('cat'))
                            map['cat'] = arr[3];
                        else if (arr[0].startsWith('timeout'))
                            map['timeout'] = parseInt(arr[4]);
                    });
                    const test = new Test((_h = map.cat) !== null && _h !== void 0 ? _h : autoTestsCatName, tests.length === 1 ? f.name : `${f.name} ${i + 1}`, () => __awaiter(this, void 0, void 0, function* () {
                        const res = yield grok.functions.eval(addNamespace(tests[i], f));
                        if (map.wait)
                            yield delay(map.wait);
                        // eslint-disable-next-line no-throw-literal
                        if (typeof res === 'boolean' && !res)
                            throw `Failed: ${tests[i]}, expected true, got ${res}`;
                    }), { skipReason: map.skip, timeout: DG.Test.isInBenchmark ? (_j = map.benchmarkTimeout) !== null && _j !== void 0 ? _j : BENCHMARK_TIMEOUT : (_k = map.timeout) !== null && _k !== void 0 ? _k : STANDART_TIMEOUT });
                    if (map.cat) {
                        const cat = map.cat;
                        if (moduleTests[cat] === undefined)
                            moduleTests[cat] = { tests: [], clear: true };
                        // only before/after can be defined in ts files tests under the category
                        if (!moduleTests[cat].tests)
                            moduleTests[cat].tests = [];
                        moduleTests[cat].tests.push(test);
                    }
                    else
                        moduleAutoTests.push(test);
                }
            }
            if (demo) {
                const wait = f.options['demoWait'] ? parseInt(f.options['demoWait']) : undefined;
                const test = new Test(demoCatName, f.friendlyName, () => __awaiter(this, void 0, void 0, function* () {
                    yield delay(300);
                    grok.shell.clearLastError();
                    yield f.apply();
                    yield delay(wait ? wait : 2000);
                    const unhandled = yield grok.shell.lastError;
                    if (unhandled)
                        throw new Error(unhandled);
                }), { skipReason: f.options['demoSkip'] });
                moduleDemo.push(test);
            }
            if (f.hasTag('semTypeDetector')) {
                let detectorsTestData = testData;
                if (f.options['testData']) {
                    detectorsTestData = yield grok.data.files.openTable(`System:AppData/${package_.nqName}/${f.options['testData']}`);
                }
                const test = new Test(detectorsCatName, f.friendlyName, () => __awaiter(this, void 0, void 0, function* () {
                    const arr = [];
                    console.log(`System:AppData/${package_.nqName}/${f.options['testData']}`);
                    for (const col of detectorsTestData.clone().columns) {
                        const res = yield f.apply([col]);
                        arr.push(res || col.semType);
                    }
                    const resArr = arr.filter((i) => i);
                    expect(resArr.length, 1);
                    if (f.options['testDataColumnName'])
                        expect(resArr[0], f.options['testDataColumnName']);
                }), { skipReason: f.options['skipTest'] });
                moduleDetectors.push(test);
            }
        }
        wasRegistered[packageId] = true;
        if (moduleAutoTests.length > 0)
            moduleTests[autoTestsCatName] = { tests: moduleAutoTests, clear: true };
        if (moduleDemo.length > 0)
            moduleTests[demoCatName] = { tests: moduleDemo, clear: true };
        if (moduleDetectors.length > 0)
            moduleTests[detectorsCatName] = { tests: moduleDetectors, clear: false };
    });
}
function redefineConsole() {
    const logs = [];
    console.log = (...args) => {
        logs.push(...args);
        stdLog(...args);
    };
    console.info = (...args) => {
        logs.push(...args);
        stdInfo(...args);
    };
    console.warn = (...args) => {
        logs.push(...args);
        stdWarn(...args);
    };
    console.error = (...args) => {
        logs.push(...args);
        stdError(...args);
    };
    return logs;
}
function resetConsole() {
    console.log = stdLog;
    console.info = stdInfo;
    console.warn = stdWarn;
    console.error = stdError;
}
export function runTests(options) {
    var _a, _b;
    var _c;
    return __awaiter(this, void 0, void 0, function* () {
        const package_ = (options === null || options === void 0 ? void 0 : options.nodeOptions) ? options.nodeOptions.package : grok.functions.getCurrentCall().func.package;
        if (!package_)
            throw new Error('Can\'t run tests outside of the package');
        const match = (_a = package_.packageOwner) === null || _a === void 0 ? void 0 : _a.match(/<([^>]*)>/);
        const packageOwner = match ? match[1] : '';
        if (package_ != undefined)
            yield initAutoTests(package_);
        const results = [];
        console.log(`Running tests`);
        console.log(options);
        options !== null && options !== void 0 ? options : (options = {});
        (_b = (_c = options).testContext) !== null && _b !== void 0 ? _b : (_c.testContext = new TestContext());
        grok.shell.clearLastError();
        const logs = redefineConsole();
        yield invokeTests(tests, options);
        for (let r of results) {
            r.result = r.result.toString().replace(/"/g, '\'');
            if (r.logs != undefined)
                r.logs = r.logs.toString().replace(/"/g, '\'');
        }
        return results;
        function invokeCategoryMethod(method, category) {
            return __awaiter(this, void 0, void 0, function* () {
                let invokationResult = undefined;
                try {
                    if (method !== undefined) {
                        yield timeout(() => __awaiter(this, void 0, void 0, function* () {
                            yield method();
                        }), 100000, `before ${category}: timeout error`);
                    }
                }
                catch (x) {
                    invokationResult = yield getResult(x);
                }
                return invokationResult;
            });
        }
        function invokeTestsInCategory(category, options, isTargetCategory) {
            var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l, _m, _o, _p, _q, _r, _s, _t;
            return __awaiter(this, void 0, void 0, function* () {
                let t = (_a = category.tests) !== null && _a !== void 0 ? _a : [];
                const res = [];
                // let memoryUsageBefore = (window?.performance as any)?.memory?.usedJSHeapSize;
                const widgetsBefore = getWidgetsCountSafe();
                if (category.clear) {
                    let skippingTests = isTargetCategory && options.skipToTest != undefined;
                    for (let i = 0; i < t.length; i++) {
                        if (t[i].options) {
                            if (((_b = t[i].options) === null || _b === void 0 ? void 0 : _b.benchmark) === undefined) {
                                if (!t[i].options)
                                    t[i].options = {};
                                t[i].options.benchmark = (_c = category.benchmarks) !== null && _c !== void 0 ? _c : false;
                            }
                        }
                        let test = t[i];
                        if (options.test)
                            if (options.test.toLowerCase() !== test.name.toLowerCase())
                                continue;
                        if (skippingTests) {
                            if ((options === null || options === void 0 ? void 0 : options.skipToTest) != undefined && test.name.toLowerCase().trim() === (options === null || options === void 0 ? void 0 : options.skipToTest.toLowerCase().trim())) {
                                // Found the target test, stop skipping after this one
                                skippingTests = false;
                            }
                            else
                                continue;
                        }
                        if (test === null || test === void 0 ? void 0 : test.options) {
                            test.options.owner = (_g = (_f = (_e = (_d = t[i].options) === null || _d === void 0 ? void 0 : _d.owner) !== null && _e !== void 0 ? _e : category === null || category === void 0 ? void 0 : category.owner) !== null && _f !== void 0 ? _f : packageOwner) !== null && _g !== void 0 ? _g : '';
                        }
                        // let isGBEnable = (window as any).gc && test.options?.skipReason == undefined;
                        // console.log(`********${isGBEnable}`);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        // memoryUsageBefore = (window?.performance as any)?.memory?.usedJSHeapSize;
                        let testRun = yield execTest(test, options === null || options === void 0 ? void 0 : options.test, logs, DG.Test.isInBenchmark ? (_j = (_h = t[i].options) === null || _h === void 0 ? void 0 : _h.benchmarkTimeout) !== null && _j !== void 0 ? _j : BENCHMARK_TIMEOUT : (_l = (_k = t[i].options) === null || _k === void 0 ? void 0 : _k.timeout) !== null && _l !== void 0 ? _l : STANDART_TIMEOUT, package_.name, options.verbose);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        if (testRun) {
                            res.push(Object.assign(Object.assign({}, testRun), { widgetsDifference: getWidgetsCountSafe() - widgetsBefore }));
                            // Return early if returnOnFail is set and test failed (but ignore failure for the skipToTest test itself)
                            if (options.returnOnFail && options.skipToTest !== test.name && !testRun.success && !testRun.skipped)
                                return res;
                        }
                        // res.push({ ...testRun, memoryDelta: (window?.performance as any)?.memory?.usedJSHeapSize - memoryUsageBefore, widgetsDelta: getWidgetsCountSafe() - widgetsBefore });
                        if (!options.nodeOptions) {
                            grok.shell.closeAll();
                            DG.Balloon.closeAll();
                        }
                    }
                }
                else {
                    let skippingTests = isTargetCategory && options.skipToTest != undefined;
                    for (let i = 0; i < t.length; i++) {
                        let test = t[i];
                        if (options.test)
                            if (options.test.toLowerCase() !== test.name.toLowerCase())
                                continue;
                        if (skippingTests) {
                            if ((options === null || options === void 0 ? void 0 : options.skipToTest) != undefined && test.name.toLowerCase().trim() === (options === null || options === void 0 ? void 0 : options.skipToTest.toLowerCase().trim())) {
                                // Found the target test, stop skipping after this one
                                skippingTests = false;
                            }
                            continue; // Skip this test (including the target)
                        }
                        if (test === null || test === void 0 ? void 0 : test.options) {
                            test.options.owner = (_q = (_p = (_o = (_m = t[i].options) === null || _m === void 0 ? void 0 : _m.owner) !== null && _o !== void 0 ? _o : category === null || category === void 0 ? void 0 : category.owner) !== null && _p !== void 0 ? _p : packageOwner) !== null && _q !== void 0 ? _q : '';
                        }
                        // let isGBEnable = (window as any).gc && test.options?.skipReason == undefined;
                        // console.log(`********${isGBEnable}`);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        // memoryUsageBefore = (window?.performance as any)?.memory?.usedJSHeapSize;
                        let testRun = yield execTest(test, options === null || options === void 0 ? void 0 : options.test, logs, DG.Test.isInBenchmark ? (_s = (_r = t[i].options) === null || _r === void 0 ? void 0 : _r.benchmarkTimeout) !== null && _s !== void 0 ? _s : BENCHMARK_TIMEOUT : (_t = t[i].options) === null || _t === void 0 ? void 0 : _t.timeout, package_.name, options.verbose);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        if (testRun) {
                            res.push(Object.assign(Object.assign({}, testRun), { widgetsDifference: getWidgetsCountSafe() - widgetsBefore }));
                            // Return early if returnOnFail is set and test failed (but ignore failure for the skipToTest test itself)
                            if (options.returnOnFail && options.skipToTest !== test.name && !testRun.success && !testRun.skipped)
                                return res;
                        }
                        // res.push({ ...testRun, memoryDelta: (window?.performance as any)?.memory?.usedJSHeapSize - memoryUsageBefore, widgetsDifference: getWidgetsCountSafe() - widgetsBefore });
                    }
                }
                return res;
            });
        }
        function getWidgetsCountSafe() {
            var _a;
            if (typeof process !== 'undefined')
                return 0;
            let length = -1;
            try {
                length = DG.Widget.getAll().length;
            }
            catch (e) {
                console.warn((_a = e.message) !== null && _a !== void 0 ? _a : e);
            }
            return length;
        }
        function invokeTests(categoriesToInvoke, options) {
            var _a, _b, _c, _d, _e;
            return __awaiter(this, void 0, void 0, function* () {
                try {
                    let skippingCategories = (options === null || options === void 0 ? void 0 : options.skipToCategory) != undefined;
                    let isTargetCategory = false;
                    for (const [key, value] of Object.entries(categoriesToInvoke)) {
                        if ((_a = options.exclude) === null || _a === void 0 ? void 0 : _a.some((c) => key.startsWith(c)))
                            continue;
                        if ((options === null || options === void 0 ? void 0 : options.category) != null && !key.toLowerCase().startsWith(`${options === null || options === void 0 ? void 0 : options.category.toLowerCase().trim()} :`) &&
                            key.toLowerCase().trim() !== (options === null || options === void 0 ? void 0 : options.category.toLowerCase().trim()))
                            continue;
                        if (skippingCategories) {
                            if (isTargetCategory)
                                skippingCategories = false;
                            else {
                                if ((options === null || options === void 0 ? void 0 : options.skipToCategory) != null && key.toLowerCase().trim() === (options === null || options === void 0 ? void 0 : options.skipToCategory.toLowerCase().trim())) {
                                    isTargetCategory = true;
                                }
                                else {
                                    // Haven't found the target category yet, keep skipping
                                    continue;
                                }
                            }
                        }
                        stdLog(`Package testing: Started {{${key}}}`);
                        //@ts-ignore
                        const skipped = (_b = value.tests) === null || _b === void 0 ? void 0 : _b.every((t) => { var _a; return (_a = t.options) === null || _a === void 0 ? void 0 : _a.skipReason; });
                        if (!skipped)
                            value.beforeStatus = yield invokeCategoryMethod(value.before, key);
                        let t = (_c = value.tests) !== null && _c !== void 0 ? _c : [];
                        if (options.stressTest) {
                            t = t.filter((e) => { var _a; return (_a = e.options) === null || _a === void 0 ? void 0 : _a.stressTest; });
                            t = shuffle(t);
                        }
                        if (((_e = (_d = options.tags) === null || _d === void 0 ? void 0 : _d.length) !== null && _e !== void 0 ? _e : 0) > 0) {
                            t = t.filter((e) => { var _a, _b; return (_b = (_a = e.options) === null || _a === void 0 ? void 0 : _a.tags) === null || _b === void 0 ? void 0 : _b.some(tag => { var _a; return ((_a = options === null || options === void 0 ? void 0 : options.tags) !== null && _a !== void 0 ? _a : []).includes(tag); }); });
                        }
                        let res;
                        if (value.beforeStatus) {
                            res = Array.from(t.map((testElem) => {
                                return {
                                    date: new Date().toISOString(),
                                    category: key,
                                    name: testElem.name,
                                    success: false,
                                    result: 'before() failed',
                                    ms: 0,
                                    skipped: false,
                                    logs: '',
                                    owner: packageOwner,
                                    package: package_.name,
                                    widgetsDifference: 0,
                                    flaking: DG.Test.isReproducing
                                };
                            }));
                            res.forEach((test) => __awaiter(this, void 0, void 0, function* () { return yield grok.shell.reportTest('package', test); }));
                        }
                        else
                            res = yield invokeTestsInCategory(value, options, skippingCategories);
                        const data = res.filter((d) => d.result != 'skipped');
                        if (!skipped)
                            value.afterStatus = yield invokeCategoryMethod(value.after, key);
                        // Clear after category
                        // grok.shell.closeAll();
                        // DG.Balloon.closeAll();
                        if (value.afterStatus) {
                            stdLog(`Package testing: Category after() {{${key}}} failed`);
                            stdLog(`Package testing: Result for {{${key}}} after: ${value.afterStatus}`);
                            data.push({
                                date: new Date().toISOString(),
                                category: key,
                                name: 'after',
                                success: false,
                                result: value.afterStatus,
                                ms: 0,
                                skipped: false,
                                logs: '',
                                owner: packageOwner,
                                package: package_.name,
                                widgetsDifference: 0,
                                flaking: DG.Test.isReproducing
                            });
                        }
                        if (value.beforeStatus) {
                            stdLog(`Package testing: Category before() {{${key}}} failed`);
                            stdLog(`Package testing: Result for {{${key}}} before: ${value.beforeStatus}`);
                            data.push({
                                date: new Date().toISOString(),
                                category: key,
                                name: 'before',
                                success: false,
                                result: value.beforeStatus,
                                ms: 0,
                                skipped: false,
                                logs: '',
                                owner: packageOwner,
                                package: package_.name,
                                widgetsDifference: 0,
                                flaking: DG.Test.isReproducing
                            });
                        }
                        results.push(...data);
                        // If returnOnFail is set and a test failed (other than skipToTest), stop processing more categories
                        if (options.returnOnFail && data.some((d) => !d.success && !d.skipped && d.name !== options.skipToTest))
                            break;
                    }
                }
                finally {
                    resetConsole();
                }
                if (options.testContext.catchUnhandled && (!DG.Test.isInBenchmark)) {
                    yield delay(1000);
                    const error = yield grok.shell.lastError;
                    if (error != undefined) {
                        const params = {
                            logs: '',
                            date: new Date().toISOString(),
                            category: 'Unhandled exceptions',
                            name: 'Exception',
                            result: error !== null && error !== void 0 ? error : '',
                            success: !error,
                            ms: 0,
                            skipped: false,
                            owner: packageOwner !== null && packageOwner !== void 0 ? packageOwner : '',
                            'package': package_.name,
                            widgetsDifference: 0
                        };
                        stdLog(`Package testing: Unhandled Exception: ${error}`);
                        results.push(Object.assign(Object.assign({}, params), { 'flaking': DG.Test.isReproducing && !error }));
                        params.package = package_.name;
                        yield grok.shell.reportTest('package', params);
                    }
                }
            });
        }
    });
}
function getResult(x) {
    return __awaiter(this, void 0, void 0, function* () {
        return `${x.toString()}\n${x.stack ? (yield DG.Logger.translateStackTrace(x.stack)) : ''}`;
    });
}
function execTest(t, predicate, logs, testTimeout, packageName, verbose) {
    var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l, _m, _o, _p;
    return __awaiter(this, void 0, void 0, function* () {
        logs.length = 0;
        let r;
        let type = 'package';
        const filter = predicate != undefined && (t.name.toLowerCase() !== predicate.toLowerCase());
        let skip = ((_a = t.options) === null || _a === void 0 ? void 0 : _a.skipReason) || filter;
        let skipReason = filter ? 'skipped' : (_b = t.options) === null || _b === void 0 ? void 0 : _b.skipReason;
        if (DG.Test.isInBenchmark && !((_c = t.options) === null || _c === void 0 ? void 0 : _c.benchmark)) {
            stdLog(`Package testing: Skipped {{${t.category}}} {{${t.name}}} doesnt available in benchmark mode`);
            return undefined;
        }
        if (!skip)
            stdLog(`Package testing: Started {{${t.category}}} {{${t.name}}}`);
        const start = Date.now();
        const startDate = new Date(start).toISOString();
        try {
            if (skip)
                r = { name: t.name, owner: (_e = (_d = t.options) === null || _d === void 0 ? void 0 : _d.owner) !== null && _e !== void 0 ? _e : '', category: t.category, logs: '', date: startDate, success: true, result: skipReason, ms: 0, skipped: true, package: packageName !== null && packageName !== void 0 ? packageName : '', flaking: DG.Test.isReproducing };
            else {
                let timeout_ = testTimeout !== null && testTimeout !== void 0 ? testTimeout : STANDART_TIMEOUT;
                if (DG.Test.isProfiling)
                    console.profile(`${t.category}: ${t.name}`);
                r = { name: t.name, owner: (_g = (_f = t.options) === null || _f === void 0 ? void 0 : _f.owner) !== null && _g !== void 0 ? _g : '', category: t.category, logs: '', date: startDate, success: true, result: (_h = (yield timeout(t.test, timeout_)).toString()) !== null && _h !== void 0 ? _h : 'OK', ms: 0, skipped: false, package: packageName !== null && packageName !== void 0 ? packageName : '', flaking: DG.Test.isReproducing };
                if (DG.Test.isProfiling) {
                    console.profileEnd(`${t.category}: ${t.name}`);
                    grok.shell.info(`Profiling of ${t.category}: ${t.name} finished \n Please ensure that you have opened DevTools (F12) / Performance panel before test starts.`);
                }
            }
        }
        catch (x) {
            stdError(x);
            r = { name: t.name, owner: (_k = (_j = t.options) === null || _j === void 0 ? void 0 : _j.owner) !== null && _k !== void 0 ? _k : '', category: t.category, logs: '', date: startDate, success: false, result: yield getResult(x), ms: 0, skipped: false, package: packageName !== null && packageName !== void 0 ? packageName : '', flaking: false };
        }
        if (((_l = t.options) === null || _l === void 0 ? void 0 : _l.isAggregated) && r.result.constructor === DG.DataFrame) {
            const col = r.result.col('success');
            if (col)
                r.success = col.stats.sum === col.length;
            if (!verbose) {
                const df = r.result;
                df.columns.remove('stack');
                df.rows.removeWhere((r) => r.get('success'));
                r.result = df;
            }
            r.result = r.result.toCsv();
        }
        r.logs = logs.join('\n');
        r.ms = Date.now() - start;
        if (!skip)
            stdLog(`Package testing: Finished {{${t.category}}} {{${t.name}}} with {{${r.success ? 'success' : 'error'}}} for ${r.ms} ms`);
        if (!r.success) {
            stdLog(`Package testing: Result for {{${t.category}}} {{${t.name}}}: ${r.result}`);
        }
        r.category = t.category;
        r.name = t.name;
        r.owner = (_o = (_m = t.options) === null || _m === void 0 ? void 0 : _m.owner) !== null && _o !== void 0 ? _o : '';
        if (!filter) {
            let params = {
                'success': r.success, 'result': r.result, 'ms': r.ms, 'date': r.date,
                'skipped': r.skipped, 'category': t.category, 'name': t.name, 'logs': r.logs, 'owner': r.owner,
                'flaking': DG.Test.isReproducing && r.success,
                'package': r.package
            };
            if (r.result.constructor == Object) {
                const res = Object.keys(r.result).reduce((acc, k) => (Object.assign(Object.assign({}, acc), { ['result.' + k]: r.result[k] })), {});
                params = Object.assign(Object.assign({}, params), res);
            }
            if (params.result instanceof DG.DataFrame)
                params.result = JSON.stringify((_p = params.result) === null || _p === void 0 ? void 0 : _p.toJson()) || '';
            yield grok.shell.reportTest(type, params);
        }
        return r;
    });
}
export function shuffle(array) {
    const newArr = array.slice();
    newArr.sort(() => Math.random() - 0.5);
    return newArr;
}
/* Waits [ms] milliseconds */
export function delay(ms) {
    return __awaiter(this, void 0, void 0, function* () {
        yield new Promise((r) => setTimeout(r, ms));
    });
}
export function awaitCheck(checkHandler, error = 'Timeout exceeded', wait = 500, interval = 50) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => {
            setTimeout(() => {
                clearInterval(intervalId);
                reject(new Error(error));
            }, wait);
            // @ts-ignore
            const intervalId = setInterval(() => {
                if (checkHandler()) {
                    clearInterval(intervalId);
                    resolve(null);
                }
            }, interval);
        });
    });
}
// Returns test execution result or an error in case of timeout
export function timeout(func, testTimeout, timeoutReason = 'EXECUTION TIMEOUT') {
    return __awaiter(this, void 0, void 0, function* () {
        let timeout = null;
        const timeoutPromise = new Promise((_, reject) => {
            timeout = setTimeout(() => {
                // eslint-disable-next-line prefer-promise-reject-errors
                reject(timeoutReason);
            }, testTimeout);
        });
        try {
            return yield Promise.race([func(), timeoutPromise]);
        }
        finally {
            if (timeout)
                clearTimeout(timeout);
        }
    });
}
export function isDialogPresent(dialogTitle) {
    const dialogs = DG.Dialog.getOpenDialogs();
    for (let i = 0; i < dialogs.length; i++) {
        if (dialogs[i].title == dialogTitle)
            return true;
    }
    return false;
}
/** Expects an asynchronous {@link action} to throw an exception. Use {@link check} to perform
 * deeper inspection of the exception if necessary.
 * @param  {function(): Promise<void>} action
 * @param  {function(any): boolean} check
 * @return {Promise<void>}
 */
export function expectExceptionAsync(action, check) {
    return __awaiter(this, void 0, void 0, function* () {
        let caught = false;
        let checked = false;
        try {
            yield action();
        }
        catch (e) {
            caught = true;
            checked = !check || check(e);
        }
        finally {
            if (!caught)
                throw new Error('An exception is expected but not thrown');
            if (!checked)
                throw new Error('An expected exception is thrown, but it does not satisfy the condition');
        }
    });
}
const catDF = DG.DataFrame.fromColumns([DG.Column.fromStrings('col', ['val1', 'val2', 'val3'])]);
/**
 * Universal test for viewers. It search viewers in DOM by tags: canvas, svg, img, input, h1, a
 * @param  {string} v Viewer name
 * @param  {_DG.DataFrame} df Dataframe to use. Should have at least 3 rows
 * @param  {boolean} options.detectSemanticTypes Specify whether to detect semantic types or not
 * @param  {boolean} options.readOnly If set to true, the dataframe will not be modified during the test
 * @param  {boolean} options.arbitraryDfTest If set to false, test on arbitrary dataframe
 * (one categorical column) will not be performed
 * @param  {object} options List of options (optional)
 * @return {Promise<void>} The test is considered successful if it completes without errors
 */
export function testViewer(v, df, options) {
    var _a;
    return __awaiter(this, void 0, void 0, function* () {
        const packageName = (_a = options === null || options === void 0 ? void 0 : options.packageName) !== null && _a !== void 0 ? _a : '';
        if (options === null || options === void 0 ? void 0 : options.detectSemanticTypes)
            yield grok.data.detectSemanticTypes(df);
        const tv = grok.shell.addTableView(df);
        try {
            //1. Open, do nothing and close
            yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded);
            //in case viewer with async rendering - wait for render to complete
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, undefined, options.awaitViewer);
            //2. Open viewer, run selection, filter, etc. and close
            if (!(options === null || options === void 0 ? void 0 : options.readOnly)) {
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, selectFilterChangeCurrent);
                if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                    yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, selectFilterChangeCurrent, options.awaitViewer);
            }
            //2. Open viewer, change options, save layout and close
            let propsAndLayout = null;
            propsAndLayout = yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, changeOptionsSaveLayout);
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                propsAndLayout = yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, changeOptionsSaveLayout, options.awaitViewer);
            //3. Load layout
            yield testViewerInternal(tv, v, packageName, grok.events.onViewLayoutApplied, loadLayout, undefined, propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.layout, { savedProps: propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.savedProps });
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                yield testViewerInternal(tv, v, packageName, grok.events.onViewLayoutApplied, loadLayout, options.awaitViewer, propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.layout, { savedProps: propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.savedProps });
            //4. Open viewer on arbitary dataset
            if ((options === null || options === void 0 ? void 0 : options.arbitraryDfTest) !== false) {
                tv.dataFrame = catDF;
                yield delay(50);
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded);
                if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                    yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, undefined, options.awaitViewer);
            }
            //5. Call postponed filtering
            yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, filterAsync);
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, filterAsync, options.awaitViewer);
        }
        finally {
            // closeAll() is handling by common test workflow
            // grok.shell.closeAll();
            // DG.Balloon.closeAll();
        }
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzIjpbInRlc3QudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7Ozs7Ozs7O0FBS0EsT0FBTyxFQUFFLFFBQVEsRUFBRSxNQUFNLG1CQUFtQixDQUFDO0FBRTdDLE9BQU8sRUFBRSx1QkFBdUIsRUFBRSxXQUFXLEVBQUUsVUFBVSxFQUFFLHlCQUF5QixFQUFFLGtCQUFrQixFQUFFLE1BQU0scUJBQXFCLENBQUM7QUFFdEksTUFBTSxnQkFBZ0IsR0FBRyxLQUFLLENBQUM7QUFDL0IsTUFBTSxpQkFBaUIsR0FBRyxRQUFRLENBQUM7QUFFbkMsTUFBTSxNQUFNLEdBQUcsT0FBTyxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDekMsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxRQUFRLEdBQUcsT0FBTyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFFN0MsTUFBTSxDQUFDLE1BQU0sS0FBSyxHQUVkLEVBQUUsQ0FBQztBQUVQLE1BQU0sZ0JBQWdCLEdBQUcsWUFBWSxDQUFDO0FBQ3RDLE1BQU0sV0FBVyxHQUFHLE1BQU0sQ0FBQztBQUMzQixNQUFNLGdCQUFnQixHQUFHLFdBQVcsQ0FBQztBQUNyQyxNQUFNLFdBQVcsR0FBRyxNQUFNLENBQUM7QUFDM0IsTUFBTSxhQUFhLEdBQStCLEVBQUUsQ0FBQztBQUNyRCxNQUFNLENBQUMsSUFBSSxlQUF1QixDQUFDO0FBRW5DLE1BQU0sS0FBVyxNQUFNLENBS3RCO0FBTEQsV0FBaUIsTUFBTTtJQUNyQixTQUFnQixPQUFPLENBQUMsS0FBVSxFQUFFLElBQWE7UUFDL0MsSUFBSSxLQUFLLElBQUksSUFBSTtZ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QUMsRUFBRSxDQUFDLENBQUM7Z0JBQ2hCLE1BQU0sa0JBQWtCLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRSxXQUFXLEVBQUUsSUFBSSxDQUFDLE1BQU0sQ0FBQyxhQUFhLENBQUMsQ0FBQztnQkFDeEUsSUFBSSxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsV0FBVztvQkFDdEIsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLGFBQWEsRUFBRSxTQUFTLEVBQUUsT0FBUSxDQUFDLFdBQVcsQ0FBQyxDQUFDO2FBQzVHO1lBRUQsNkJBQTZCO1lBQzdCLE1BQU0sa0JBQWtCLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRSxXQUFXLEVBQUUsSUFBSSxDQUFDLE1BQU0sQ0FBQyxhQUFhLEVBQUUsV0FBVyxDQUFDLENBQUM7WUFDckYsSUFBSSxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsV0FBVztnQkFDdEIsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLGFBQWEsRUFBRSxXQUFXLEVBQUUsT0FBUSxDQUFDLFdBQVcsQ0FBQyxDQUFDO1NBRTlHO2dCQUFTO1lBQ1IsaURBQWlEO1lBQ2pELHlCQUF5QjtZQUN6Qix5QkFBeUI7U0FDMUI7O0NBQ0YiLCJzb3VyY2VzQ29udGVudCI6WyJpbXBvcnQgdHlwZSAqIGFzIF9ncm9rIGZyb20gJ2RhdGFncm9rLWFwaS9ncm9rJztcclxuaW1wb3J0IHR5cGUgKiBhcyBfREcgZnJvbSAnZGF0YWdyb2stYXBpL2RnJztcclxuZGVjbGFyZSBsZXQgZ3JvazogdHlwZW9mIF9ncm9rLCBERzogdHlwZW9mIF9ERztcclxuXHJcbmltcG9ydCB7IE9ic2VydmFibGUgfSBmcm9tICdyeGpzJztcclxuaW1wb3J0IHsgdGVzdERhdGEgfSBmcm9tICcuL2RhdGFmcmFtZS11dGlscyc7XHJcbmltcG9ydCBUaW1lb3V0ID0gTm9kZUpTLlRpbWVvdXQ7XHJcbmltcG9ydCB7IGNoYW5nZU9wdGlvbnNTYXZlTGF5b3V0LCBmaWx0ZXJBc3luYywgbG9hZExheW91dCwgc2VsZWN0RmlsdGVyQ2hhbmdlQ3VycmVudCwgdGVzdFZpZXdlckludGVybmFsIH0gZnJvbSAnLi90ZXN0LXZpZXdlci11dGlscyc7XHJcblxyXG5jb25zdCBTVEFOREFSVF9USU1FT1VUID0gMzAwMDA7XHJcbmNvbnN0IEJFTkNITUFSS19USU1FT1VUID0gMTA4MDAwMDA7XHJcblxyXG5jb25zdCBzdGRMb2cgPSBjb25zb2xlLmxvZy5iaW5kKGNvbnNvbGUpO1xyXG5jb25zdCBzdGRJbmZvID0gY29uc29sZS5pbmZvLmJpbmQoY29uc29sZSk7XHJcbmNvbnN0IHN0ZFdhcm4gPSBjb25zb2xlLndhcm4uYmluZChjb25zb2xlKTtcclxuY29uc3Qgc3RkRXJyb3IgPSBjb25zb2xlLmVycm9yLmJpbmQoY29uc29sZSk7XHJcblxyXG5leHBvcnQgY29uc3QgdGVzdHM6IHtcclxuICBba2V5OiBzdHJpbmddOiBDYXRlZ29yeVxyXG59ID0ge307XHJcblxyXG5jb25zdCBhdXRvVGVzdHNDYXROYW1lID0gJ0F1dG8gVGVzdHMnO1xyXG5jb25zdCBkZW1vQ2F0TmFtZSA9ICdEZW1vJztcclxuY29uc3QgZGV0ZWN0b3JzQ2F0TmFtZSA9ICdEZXRlY3RvcnMnO1xyXG5jb25zdCBjb3JlQ2F0TmFtZSA9ICdDb3JlJztcclxuY29uc3Qgd2FzUmVnaXN0ZXJlZDogeyBba2V5OiBzdHJpbmddOiBib29sZWFuIH0gPSB7fTtcclxuZXhwb3J0IGxldCBjdXJyZW50Q2F0ZWdvcnk6IHN0cmluZztcclxuXHJcbmV4cG9ydCBuYW1lc3BhY2UgYXNzdXJlIHtcclxuICBleHBvcnQgZnVuY3Rpb24gbm90TnVsbCh2YWx1ZTogYW55LCBuYW1lPzogc3RyaW5nKSB7XHJcbiAgICBpZiAodmFsdWUgPT0gbnVsbClcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKGAke25hbWUgPT0gbnVsbCA/ICdWYWx1ZScgOiBuYW1lfSBub3QgZGVmaW5lZGApO1xyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGludGVyZmFjZSBUZXN0T3B0aW9ucyB7XHJcbiAgdGltZW91dD86IG51bWJlcjtcclxuICBiZW5jaG1hcmtXYXJuVGltZW91dD86IG51bWJlcjtcclxuICBiZW5jaG1hcmtUaW1lb3V0PzogbnVtYmVyO1xyXG4gIHVuaGFuZGxlZEV4Y2VwdGlvblRpbWVvdXQ/OiBudW1iZXI7XHJcbiAgc2tpcFJlYXNvbj86IHN0cmluZztcclxuICBpc0FnZ3JlZ2F0ZWQ/OiBib29sZWFuO1xyXG4gIGJlbmNobWFyaz86IGJvb2xlYW47XHJcbiAgc3RyZXNzVGVzdD86IGJvb2xlYW47XHJcbiAgb3duZXI/OiBzdHJpbmc7XHJcbiAgdGFncz86IHN0cmluZ1tdO1xyXG59XHJcblxyXG5leHBvcnQgaW50ZXJmYWNlIFRlc3RSZXN1bHQge1xyXG4gIGRhdGU6IHN0cmluZztcclxuICBjYXRlZ29yeTogc3RyaW5nO1xyXG4gIG5hbWU6IHN0cmluZztcclxuICBzdWNjZXNzOiBib29sZWFuO1xyXG4gIHJlc3VsdDogYW55O1xyXG4gIG1zOiBudW1iZXI7XHJcbiAgc2tpcHBlZDogYm9vbGVhbjtcclxuICBsb2dzOiBzdHJpbmc7XHJcbiAgb3duZXI6IHN0cmluZztcclxuICBwYWNrYWdlOiBzdHJpbmc7XHJcbiAgZmxha2luZzogYm9vbGVhbjtcclxufVxyXG5cclxuXHJcbmV4cG9ydCBpbnRlcmZhY2UgVGVzdFJlc3VsdEV4dGVuZGVkIGV4dGVuZHMgVGVzdFJlc3VsdHtcclxuICB3aWRnZXRzRGlmZmVyZW5jZTogbnVtYmVyO1xyXG59XHJcblxyXG5leHBvcnQgaW50ZXJmYWNlIENhdGVnb3J5T3B0aW9ucyB7XHJcbiAgY2xlYXI/OiBib29sZWFuO1xyXG4gIHRpbWVvdXQ/OiBudW1iZXI7XHJcbiAgYmVuY2htYXJrcz86IGJvb2xlYW47XHJcbiAgc3RyZXNzVGVzdHM/OiBib29sZWFuO1xyXG4gIG93bmVyPzogc3RyaW5nO1xyXG59XHJcblxyXG5leHBvcnQgY2xhc3MgVGVzdENvbnRleHQge1xyXG4gIHN0cmVzc1Rlc3Q/OiBib29sZWFuO1xyXG4gIGNhdGNoVW5oYW5kbGVkID0gdHJ1ZTtcclxuICByZXBvcnQgPSBmYWxzZTtcclxuICByZXR1cm5PbkZhaWwgPSBmYWxzZTtcclxuXHJcbiAgY29uc3RydWN0b3IoY2F0Y2hVbmhhbmRsZWQ/OiBib29sZWFuLCByZXBvcnQ/OiBib29sZWFuLCByZXR1cm5PbkZhaWw/OiBib29sZWFuKSB7XHJcbiAgICBpZiAoY2F0Y2hVbmhhbmRsZWQgIT09IHVuZGVmaW5lZCkgdGhpcy5jYXRjaFVuaGFuZGxlZCA9IGNhdGNoVW5oYW5kbGVkO1xyXG4gICAgaWYgKHJlcG9ydCAhPT0gdW5kZWZpbmVkKSB0aGlzLnJlcG9ydCA9IHJlcG9ydDtcclxuICAgIGlmIChyZXR1cm5PbkZhaWwgIT09IHVuZGVmaW5lZCkgdGhpcy5yZXR1cm5PbkZhaWwgPSByZXR1cm5PbkZhaWw7XHJcbiAgfTtcclxufVxyXG5cclxuZXhwb3J0IGNsYXNzIFRlc3Qge1xyXG4gIHRlc3Q6ICgpID0+IFByb21pc2U8YW55PjtcclxuICBuYW1lOiBzdHJpbmc7XHJcbiAgY2F0ZWdvcnk6IHN0cmluZztcclxuICBvcHRpb25zPzogVGVzdE9wdGlvbnM7XHJcblxyXG4gIGNvbnN0cnVjdG9yKGNhdGVnb3J5OiBzdHJpbmcsIG5hbWU6IHN0cmluZywgdGVzdDogKCkgPT4gUHJvbWlzZTxhbnk+LCBvcHRpb25zPzogVGVzdE9wdGlvbnMpIHtcclxuICAgIHRoaXMuY2F0ZWdvcnkgPSBjYXRlZ29yeTtcclxuICAgIHRoaXMubmFtZSA9IG5hbWU7XHJcbiAgICBvcHRpb25zID8/PSB7fTtcclxuICAgIG9wdGlvbnMudGltZW91dCA/Pz0gU1RBTkRBUlRfVElNRU9VVDtcclxuICAgIHRoaXMub3B0aW9ucyA9IG9wdGlvbnM7XHJcbiAgICB0aGlzLnRlc3QgPSBhc3luYyAoKTogUHJvbWlzZTxhbnk+ID0+IHtcclxuICAgICAgcmV0dXJuIG5ldyBQcm9taXNlKGFzeW5jIChyZXNvbHZlLCByZWplY3QpID0+IHtcclxuICAgICAgICBsZXQgcmVzdWx0ID0gJyc7XHJcbiAgICAgICAgdHJ5IHtcclxuICAgICAgICAgIGlmIChERy5UZXN0LmlzSW5EZWJ1ZylcclxuICAgICAgICAgICAgZGVidWdnZXI7XHJcblxyXG4gICAgICAgICAgbGV0IHJlcyA9IGF3YWl0IHRlc3QoKTtcclxuICAgICAgICAgIHRyeSB7XHJcbiAgICAgICAgICAgIHJlc3VsdCA9IHJlcz8udG9TdHJpbmcoKSA/PyAnJztcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGNhdGNoIChlKSB7XHJcbiAgICAgICAgICAgIHJlc3VsdCA9ICdDYW5cXCd0IGNvbnZlcnQgdGVzdFxcJ3MgcmVzdWx0IHRvIHN0cmluZyc7XHJcbiAgICAgICAgICAgIGNvbnNvbGUuZXJyb3IoYENhblxcJ3QgY29udmVydCB0ZXN0XFwncyByZXN1bHQgdG8gc3RyaW5nIGluIHRoZSAke3RoaXMuY2F0ZWdvcnl9OiR7dGhpcy5uYW1lfSB0ZXN0YCk7XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgfSBjYXRjaCAoZTogYW55KSB7XHJcbiAgICAgICAgICByZWplY3QoZSk7XHJcbiAgICAgICAgfVxyXG4gICAgICAgIHJlc29sdmUocmVzdWx0KTtcclxuICAgICAgfSk7XHJcbiAgICB9O1xyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGNsYXNzIENhdGVnb3J5IHtcclxuICB0ZXN0cz86IFRlc3RbXTtcclxuICBiZWZvcmU/OiAoKSA9PiBQcm9taXNlPHZvaWQ+O1xyXG4gIGFmdGVyPzogKCkgPT4gUHJvbWlzZTx2b2lkPjtcclxuXHJcbiAgYmVmb3JlU3RhdHVzPzogc3RyaW5nO1xyXG4gIGFmdGVyU3RhdHVzPzogc3RyaW5nO1xyXG4gIGNsZWFyPzogYm9vbGVhbjtcclxuICB0aW1lb3V0PzogbnVtYmVyO1xyXG4gIGJlbmNobWFya3M/OiBib29sZWFuO1xyXG4gIGJlbmNobWFya1RpbWVvdXQ/OiBudW1iZXI7XHJcbiAgc3RyZXNzVGVzdHM/OiBib29sZWFuO1xyXG4gIG93bmVyPzogc3RyaW5nO1xyXG59XHJcblxyXG5leHBvcnQgY2xhc3MgTm9kZVRlc3RFeGVjdXRpb25PcHRpb25zIHtcclxuICBwYWNrYWdlITogX0RHLlBhY2thZ2U7XHJcbn1cclxuXHJcbmV4cG9ydCBjbGFzcyBUZXN0RXhlY3V0aW9uT3B0aW9ucyB7XHJcbiAgY2F0ZWdvcnk/OiBzdHJpbmc7XHJcbiAgdGVzdD86IHN0cmluZztcclxuICB0ZXN0Q29udGV4dD86IFRlc3RDb250ZXh0O1xyXG4gIGV4Y2x1ZGU/OiBzdHJpbmdbXTtcclxuICB2ZXJib3NlPzogYm9vbGVhbjtcclxuICBzdHJlc3NUZXN0PzogYm9vbGVhbjtcclxuICB0YWdzPzogc3RyaW5nW107XHJcbiAgbm9kZU9wdGlvbnM/OiBOb2RlVGVzdEV4ZWN1dGlvbk9wdGlvbnM7XHJcbiAgc2tpcFRvQ2F0ZWdvcnk/OiBzdHJpbmc7XHJcbiAgc2tpcFRvVGVzdD86IHN0cmluZztcclxuICByZXR1cm5PbkZhaWw/OiBib29sZWFuO1xyXG59XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdGVzdEV2ZW50PFQ+KGV2ZW50OiBPYnNlcnZhYmxlPFQ+LFxyXG4gIGhhbmRsZXI6IChhcmdzOiBUKSA9PiB2b2lkLCB0cmlnZ2VyOiAoKSA9PiB2b2lkLCBtczogbnVtYmVyID0gMCwgcmVhc29uOiBzdHJpbmcgPSBgdGltZW91dGBcclxuKTogUHJvbWlzZTxhbnk+IHtcclxuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUsIHJlamVjdCkgPT4ge1xyXG4gICAgY29uc3Qgc3ViID0gZXZlbnQuc3Vic2NyaWJlKChhcmdzOiBUKSA9PiB7XHJcbiAgICAgIHRyeSB7XHJcbiAgICAgICAgaGFuZGxlcihhcmdzKTtcclxuICAgICAgICByZXNvbHZlKCdPSycpO1xyXG4gICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgcmVqZWN0KGUpO1xyXG4gICAgICB9IGZpbmFsbHkge1xyXG4gICAgICAgIHN1Yi51bnN1YnNjcmliZSgpO1xyXG4gICAgICAgIGNsZWFyVGltZW91dCh0aW1lb3V0KTtcclxuICAgICAgfVxyXG4gICAgfSk7XHJcbiAgICBjb25zdCB0aW1lb3V0ID0gc2V0VGltZW91dCgoKSA9PiB7XHJcbiAgICAgIHN1Yi51bnN1YnNjcmliZSgpO1xyXG4gICAgICAvLyBlc2xpbnQtZGlzYWJsZS1uZXh0LWxpbmUgcHJlZmVyLXByb21pc2UtcmVqZWN0LWVycm9yc1xyXG4gICAgICByZWplY3QocmVhc29uKTtcclxuICAgIH0sIG1zKTtcclxuICAgIHRyaWdnZXIoKTtcclxuICB9KTtcclxufVxyXG5cclxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRlc3RFdmVudEFzeW5jPFQ+KGV2ZW50OiBPYnNlcnZhYmxlPFQ+LFxyXG4gIGhhbmRsZXI6IChhcmdzOiBUKSA9PiBQcm9taXNlPHZvaWQ+LCB0cmlnZ2VyOiAoKSA9PiB2b2lkLCBtczogbnVtYmVyID0gMCwgcmVhc29uOiBzdHJpbmcgPSBgdGltZW91dGBcclxuKTogUHJvbWlzZTxhbnk+IHtcclxuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUsIHJlamVjdCkgPT4ge1xyXG4gICAgY29uc3Qgc3ViID0gZXZlbnQuc3Vic2NyaWJlKChhcmdzOiBUKSA9PiB7XHJcbiAgICAgIGhhbmRsZXIoYXJncykudGhlbigoKSA9PiB7XHJcbiAgICAgICAgcmVzb2x2ZSgnT0snKTtcclxuICAgICAgfSkuY2F0Y2goKGUpID0+IHtcclxuICAgICAgICByZWplY3QoZSk7XHJcbiAgICAgIH0pLmZpbmFsbHkoKCkgPT4ge1xyXG4gICAgICAgIHN1Yi51bnN1YnNjcmliZSgpO1xyXG4gICAgICAgIGNsZWFyVGltZW91dCh0aW1lb3V0KTtcclxuICAgICAgfSk7XHJcbiAgICB9KTtcclxuICAgIGNvbnN0IHRpbWVvdXQgPSBzZXRUaW1lb3V0KCgpID0+IHtcclxuICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XHJcbiAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBwcmVmZXItcHJvbWlzZS1yZWplY3QtZXJyb3JzXHJcbiAgICAgIHJlamVjdChyZWFzb24pO1xyXG4gICAgfSwgbXMpO1xyXG4gICAgdHJpZ2dlcigpO1xyXG4gIH0pO1xyXG59XHJcblxyXG5leHBvcnQgZnVuY3Rpb24gdGVzdChuYW1lOiBzdHJpbmcsIHRlc3Q6ICgpID0+IFByb21pc2U8YW55Piwgb3B0aW9ucz86IFRlc3RPcHRpb25zKTogdm9pZCB7XHJcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPT0gdW5kZWZpbmVkKVxyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XSA9IHt9O1xyXG4gIGlmICh0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLnRlc3RzID09IHVuZGVmaW5lZClcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGVzdHMgPSBbXTtcclxuICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLnRlc3RzIS5wdXNoKG5ldyBUZXN0KGN1cnJlbnRDYXRlZ29yeSwgbmFtZSwgdGVzdCwgb3B0aW9ucykpO1xyXG59XHJcblxyXG4vKiBUZXN0cyB0d28gb2JqZWN0cyBmb3IgZXF1YWxpdHksIHRocm93cyBhbiBleGNlcHRpb24gaWYgdGhleSBhcmUgbm90IGVxdWFsLiAqL1xyXG5leHBvcnQgZnVuY3Rpb24gZXhwZWN0KGFjdHVhbDogYW55LCBleHBlY3RlZDogYW55ID0gdHJ1ZSwgZXJyb3I/OiBzdHJpbmcpOiB2b2lkIHtcclxuICBpZiAoZXJyb3IpXHJcbiAgICBlcnJvciA9IGAke2Vycm9yfSwgYDtcclxuICBlbHNlIGVycm9yID0gJyc7XHJcbiAgaWYgKGFjdHVhbCAhPT0gZXhwZWN0ZWQpXHJcbiAgICB0aHJvdyBuZXcgRXJyb3IoYCR7ZXJyb3J9RXhwZWN0ZWQgXCIke2V4cGVjdGVkfVwiLCBnb3QgXCIke2FjdHVhbH1cImApO1xyXG59XHJcblxyXG5leHBvcnQgZnVuY3Rpb24gZXhwZWN0RmxvYXQoYWN0dWFsOiBudW1iZXIsIGV4cGVjdGVkOiBudW1iZXIsIHRvbGVyYW5jZSA9IDAuMDAxLCBlcnJvcj86IHN0cmluZyk6IHZvaWQge1xyXG4gIGlmICgoYWN0dWFsID09PSBOdW1iZXIuUE9TSVRJVkVfSU5GSU5JVFkgJiYgZXhwZWN0ZWQgPT09IE51bWJlci5QT1NJVElWRV9JTkZJTklUWSkgfHxcclxuICAgIChhY3R1YWwgPT09IE51bWJlci5ORUdBVElWRV9JTkZJTklUWSAmJiBleHBlY3RlZCA9PT0gTnVtYmVyLk5FR0FUSVZFX0lORklOSVRZKSB8fFxyXG4gICAgKGFjdHVhbCA9PT0gTnVtYmVyLk5hTiAmJiBleHBlY3RlZCA9PT0gTnVtYmVyLk5hTikgfHwgKGlzTmFOKGFjdHVhbCkgJiYgaXNOYU4oZXhwZWN0ZWQpKSlcclxuICAgIHJldHVybjtcclxuICBjb25zdCBhcmVFcXVhbCA9IE1hdGguYWJzKGFjdHVhbCAtIGV4cGVjdGVkKSA8IHRvbGVyYW5jZTtcclxuICBleHBlY3QoYXJlRXF1YWwsIHRydWUsIGAke2Vycm9yID8/ICcnfSAodG9sZXJhbmNlID0gJHt0b2xlcmFuY2V9KWApO1xyXG4gIGlmICghYXJlRXF1YWwpXHJcbiAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkICR7ZXhwZWN0ZWR9LCBnb3QgJHthY3R1YWx9ICh0b2xlcmFuY2UgPSAke3RvbGVyYW5jZX0pYCk7XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiBleHBlY3RUYWJsZShhY3R1YWw6IF9ERy5EYXRhRnJhbWUsIGV4cGVjdGVkOiBfREcuRGF0YUZyYW1lLCBlcnJvcj86IHN0cmluZyk6IHZvaWQge1xyXG4gIGNvbnN0IGV4cGVjdGVkUm93Q291bnQgPSBleHBlY3RlZC5yb3dDb3VudDtcclxuICBjb25zdCBhY3R1YWxSb3dDb3VudCA9IGFjdHVhbC5yb3dDb3VudDtcclxuICBleHBlY3QoYWN0dWFsUm93Q291bnQsIGV4cGVjdGVkUm93Q291bnQsIGAke2Vycm9yID8/ICcnfSwgcm93IGNvdW50YCk7XHJcblxyXG4gIGZvciAoY29uc3QgY29sdW1uIG9mIGV4cGVjdGVkLmNvbHVtbnMpIHtcclxuICAgIGNvbnN0IGFjdHVhbENvbHVtbiA9IGFjdHVhbC5jb2x1bW5zLmJ5TmFtZShjb2x1bW4ubmFtZSk7XHJcbiAgICBpZiAoYWN0dWFsQ29sdW1uID09IG51bGwpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihgQ29sdW1uICR7Y29sdW1uLm5hbWV9IG5vdCBmb3VuZGApO1xyXG4gICAgaWYgKGFjdHVhbENvbHVtbi50eXBlICE9IGNvbHVtbi50eXBlKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYENvbHVtbiAke2NvbHVtbi5uYW1lfSB0eXBlIGV4cGVjdGVkICR7Y29sdW1uLnR5cGV9IGdvdCAke2FjdHVhbENvbHVtbi50eXBlfWApO1xyXG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBleHBlY3RlZFJvd0NvdW50OyBpKyspIHtcclxuICAgICAgY29uc3QgdmFsdWUgPSBjb2x1bW4uZ2V0KGkpO1xyXG4gICAgICBjb25zdCBhY3R1YWxWYWx1ZSA9IGFjdHVhbENvbHVtbi5nZXQoaSk7XHJcbiAgICAgIGlmIChjb2x1bW4udHlwZSA9PSBERy5UWVBFLkZMT0FUKVxyXG4gICAgICAgIGV4cGVjdEZsb2F0KGFjdHVhbFZhbHVlLCB2YWx1ZSwgMC4wMDAxLCBlcnJvcik7XHJcbiAgICAgIGVsc2UgaWYgKGNvbHVtbi50eXBlID09IERHLlRZUEUuREFURV9USU1FKVxyXG4gICAgICAgIGV4cGVjdChhY3R1YWxWYWx1ZS5pc1NhbWUodmFsdWUpLCB0cnVlLCBlcnJvcik7XHJcbiAgICAgIGVsc2VcclxuICAgICAgICBleHBlY3QoYWN0dWFsVmFsdWUsIHZhbHVlLCBlcnJvcik7XHJcbiAgICB9XHJcbiAgfVxyXG59XHJcblxyXG5leHBvcnQgZnVuY3Rpb24gZXhwZWN0T2JqZWN0KGFjdHVhbDogeyBba2V5OiBzdHJpbmddOiBhbnkgfSwgZXhwZWN0ZWQ6IHsgW2tleTogc3RyaW5nXTogYW55IH0pIHtcclxuICBmb3IgKGNvbnN0IFtleHBlY3RlZEtleSwgZXhwZWN0ZWRWYWx1ZV0gb2YgT2JqZWN0LmVudHJpZXMoZXhwZWN0ZWQpKSB7XHJcbiAgICBpZiAoIWFjdHVhbC5oYXNPd25Qcm9wZXJ0eShleHBlY3RlZEtleSkpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihgRXhwZWN0ZWQgcHJvcGVydHkgXCIke2V4cGVjdGVkS2V5fVwiIG5vdCBmb3VuZGApO1xyXG5cclxuICAgIGNvbnN0IGFjdHVhbFZhbHVlID0gYWN0dWFsW2V4cGVjdGVkS2V5XTtcclxuICAgIGlmIChhY3R1YWxWYWx1ZSBpbnN0YW5jZW9mIEFycmF5ICYmIGV4cGVjdGVkVmFsdWUgaW5zdGFuY2VvZiBBcnJheSlcclxuICAgICAgZXhwZWN0QXJyYXkoYWN0dWFsVmFsdWUsIGV4cGVjdGVkVmFsdWUpO1xyXG4gICAgZWxzZSBpZiAoYWN0dWFsVmFsdWUgaW5zdGFuY2VvZiBPYmplY3QgJiYgZXhwZWN0ZWRWYWx1ZSBpbnN0YW5jZW9mIE9iamVjdClcclxuICAgICAgZXhwZWN0T2JqZWN0KGFjdHVhbFZhbHVlLCBleHBlY3RlZFZhbHVlKTtcclxuICAgIGVsc2UgaWYgKE51bWJlci5pc0Zpbml0ZShhY3R1YWxWYWx1ZSkgJiYgTnVtYmVyLmlzRmluaXRlKGV4cGVjdGVkVmFsdWUpKVxyXG4gICAgICBleHBlY3RGbG9hdChhY3R1YWxWYWx1ZSwgZXhwZWN0ZWRWYWx1ZSk7XHJcbiAgICBlbHNlIGlmIChhY3R1YWxWYWx1ZSAhPSBleHBlY3RlZFZhbHVlKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkICgke2V4cGVjdGVkVmFsdWV9KSBmb3Iga2V5ICcke2V4cGVjdGVkS2V5fScsIGdvdCAoJHthY3R1YWxWYWx1ZX0pYCk7XHJcbiAgfVxyXG59XHJcblxyXG5leHBvcnQgZnVuY3Rpb24gZXhwZWN0QXJyYXkoYWN0dWFsOiBBcnJheUxpa2U8YW55PiwgZXhwZWN0ZWQ6IEFycmF5TGlrZTxhbnk+KSB7XHJcbiAgY29uc3QgYWN0dWFsTGVuZ3RoID0gYWN0dWFsLmxlbmd0aDtcclxuICBjb25zdCBleHBlY3RlZExlbmd0aCA9IGV4cGVjdGVkLmxlbmd0aDtcclxuXHJcbiAgaWYgKGFjdHVhbExlbmd0aCAhPSBleHBlY3RlZExlbmd0aCkge1xyXG4gICAgdGhyb3cgbmV3IEVycm9yKGBBcnJheXMgYXJlIG9mIGRpZmZlcmVudCBsZW5ndGg6IGFjdHVhbCBhcnJheSBsZW5ndGggaXMgJHthY3R1YWxMZW5ndGh9IGAgK1xyXG4gICAgICBgYW5kIGV4cGVjdGVkIGFycmF5IGxlbmd0aCBpcyAke2V4cGVjdGVkTGVuZ3RofWApO1xyXG4gIH1cclxuXHJcbiAgZm9yIChsZXQgaSA9IDA7IGkgPCBhY3R1YWxMZW5ndGg7IGkrKykge1xyXG4gICAgaWYgKGFjdHVhbFtpXSBpbnN0YW5jZW9mIEFycmF5ICYmIGV4cGVjdGVkW2ldIGluc3RhbmNlb2YgQXJyYXkpXHJcbiAgICAgIGV4cGVjdEFycmF5KGFjdHVhbFtpXSwgZXhwZWN0ZWRbaV0pO1xyXG4gICAgZWxzZSBpZiAoYWN0dWFsW2ldIGluc3RhbmNlb2YgT2JqZWN0ICYmIGV4cGVjdGVkW2ldIGluc3RhbmNlb2YgT2JqZWN0KVxyXG4gICAgICBleHBlY3RPYmplY3QoYWN0dWFsW2ldLCBleHBlY3RlZFtpXSk7XHJcbiAgICBlbHNlIGlmIChhY3R1YWxbaV0gIT0gZXhwZWN0ZWRbaV0pXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihgRXhwZWN0ZWQgJHtleHBlY3RlZFtpXX0gYXQgcG9zaXRpb24gJHtpfSwgZ290ICR7YWN0dWFsW2ldfWApO1xyXG4gIH1cclxufVxyXG5cclxuLyogRGVmaW5lcyBhIHRlc3Qgc3VpdGUuICovXHJcbmV4cG9ydCBmdW5jdGlvbiBjYXRlZ29yeShjYXRlZ29yeTogc3RyaW5nLCB0ZXN0c186ICgpID0+IHZvaWQsIG9wdGlvbnM/OiBDYXRlZ29yeU9wdGlvbnMpOiB2b2lkIHtcclxuICBjdXJyZW50Q2F0ZWdvcnkgPSBjYXRlZ29yeTtcclxuICB0ZXN0c18oKTtcclxuICBpZiAodGVzdHNbY3VycmVudENhdGVnb3J5XSkge1xyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5jbGVhciA9IG9wdGlvbnM/LmNsZWFyID8/IHRydWU7XHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLnRpbWVvdXQgPSBvcHRpb25zPy50aW1lb3V0O1xyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5iZW5jaG1hcmtzID0gb3B0aW9ucz8uYmVuY2htYXJrcztcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uc3RyZXNzVGVzdHMgPSBvcHRpb25zPy5zdHJlc3NUZXN0cztcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0ub3duZXIgPSBvcHRpb25zPy5vd25lcjtcclxuICB9XHJcbn1cclxuXHJcbi8qIERlZmluZXMgYSBmdW5jdGlvbiB0byBiZSBleGVjdXRlZCBiZWZvcmUgdGhlIHRlc3RzIGluIHRoaXMgY2F0ZWdvcnkgYXJlIGV4ZWN1dGVkLiAqL1xyXG5leHBvcnQgZnVuY3Rpb24gYmVmb3JlKGJlZm9yZTogKCkgPT4gUHJvbWlzZTx2b2lkPik6IHZvaWQge1xyXG4gIGlmICh0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldID09IHVuZGVmaW5lZClcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPSB7fTtcclxuICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLmJlZm9yZSA9IGJlZm9yZTtcclxufVxyXG5cclxuLyogRGVmaW5lcyBhIGZ1bmN0aW9uIHRvIGJlIGV4ZWN1dGVkIGFmdGVyIHRoZSB0ZXN0cyBpbiB0aGlzIGNhdGVnb3J5IGFyZSBleGVjdXRlZC4gKi9cclxuZXhwb3J0IGZ1bmN0aW9uIGFmdGVyKGFmdGVyOiAoKSA9PiBQcm9taXNlPHZvaWQ+KTogdm9pZCB7XHJcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPT0gdW5kZWZpbmVkKVxyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XSA9IHt9O1xyXG4gIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uYWZ0ZXIgPSBhZnRlcjtcclxufVxyXG5cclxuZnVuY3Rpb24gYWRkTmFtZXNwYWNlKHM6IHN0cmluZywgZjogX0RHLkZ1bmMpOiBzdHJpbmcge1xyXG4gIHJldHVybiBzLnJlcGxhY2UobmV3IFJlZ0V4cChmLm5hbWUsICdnaScpLCBmLm5xTmFtZSk7XHJcbn1cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBpbml0QXV0b1Rlc3RzKHBhY2thZ2VfOiBfREcuUGFja2FnZSwgbW9kdWxlPzogYW55KSB7XHJcbiAgY29uc3QgcGFja2FnZUlkID0gcGFja2FnZV8uaWQ7XHJcbiAgaWYgKHdhc1JlZ2lzdGVyZWRbcGFja2FnZUlkXSkgcmV0dXJuO1xyXG4gIGNvbnN0IG1vZHVsZVRlc3RzID0gbW9kdWxlID8gbW9kdWxlLnRlc3RzIDogdGVzdHM7XHJcbiAgaWYgKHBhY2thZ2VfLm5hbWUgPT09ICdEZXZUb29scycgfHwgKCEhbW9kdWxlICYmIG1vZHVsZS5fcGFja2FnZS5uYW1lID09PSAnRGV2VG9vbHMnKSkge1xyXG4gICAgZm9yIChjb25zdCBmIG9mICg8YW55PndpbmRvdykuZGFydFRlc3RzKSB7XHJcbiAgICAgIGNvbnN0IGFyciA9IGYubmFtZS5zcGxpdCgvXFxzKlxcfFxccyohL2cpO1xyXG4gICAgICBsZXQgbmFtZSA9IGFyci5wb3AoKSA/PyBmLm5hbWU7XHJcbiAgICAgIGxldCBjYXQgPSBhcnIubGVuZ3RoID8gY29yZUNhdE5hbWUgKyAnOiAnICsgYXJyLmpvaW4oJzogJykgOiBjb3JlQ2F0TmFtZTtcclxuICAgICAgbGV0IGZ1bGxOYW1lOiBzdHJpbmdbXSA9IG5hbWUuc3BsaXQoJyB8ICcpO1xyXG4gICAgICBuYW1lID0gZnVsbE5hbWVbZnVsbE5hbWUubGVuZ3RoIC0gMV07XHJcbiAgICAgIGZ1bGxOYW1lLnVuc2hpZnQoY2F0KTtcclxuICAgICAgZnVsbE5hbWUucG9wKCk7XHJcbiAgICAgIGNhdCA9IGZ1bGxOYW1lLmpvaW4oJzogJyk7XHJcbiAgICAgIGlmIChtb2R1bGVUZXN0c1tjYXRdID09PSB1bmRlZmluZWQpXHJcbiAgICAgICAgbW9kdWxlVGVzdHNbY2F0XSA9IHsgdGVzdHM6IFtdLCBjbGVhcjogdHJ1ZSB9O1xyXG4gICAgICBtb2R1bGVUZXN0c1tjYXRdLnRlc3RzLnB1c2gobmV3IFRlc3QoY2F0LCBuYW1lLCBmLnRlc3QsIHsgaXNBZ2dyZWdhdGVkOiBmYWxzZSwgdGltZW91dDogZi5vcHRpb25zPy50aW1lb3V0ID8/IFNUQU5EQVJUX1RJTUVPVVQsIHNraXBSZWFzb246IGYub3B0aW9ucz8uc2tpcFJlYXNvbiwgb3duZXI6IGYub3B0aW9ucz8ub3duZXIsIGJlbmNobWFyazogZi5vcHRpb25zPy5iZW5jaG1hcmsgPz8gZmFsc2UgfSkpO1xyXG4gICAgfVxyXG4gIH1cclxuICBjb25zdCBtb2R1bGVBdXRvVGVzdHMgPSBbXTtcclxuICBjb25zdCBtb2R1bGVEZW1vID0gW107XHJcbiAgY29uc3QgbW9kdWxlRGV0ZWN0b3JzID0gW107XHJcbiAgY29uc3QgcGFja0Z1bmN0aW9ucyA9IGF3YWl0IGdyb2suZGFwaS5mdW5jdGlvbnMuZmlsdGVyKGBwYWNrYWdlLmlkID0gXCIke3BhY2thZ2VJZH1cImApLmxpc3QoKTtcclxuICBjb25zdCByZWcgPSBuZXcgUmVnRXhwKC9za2lwOlxccyooW14sXFxzXSspfHdhaXQ6XFxzKihcXGQrKXxjYXQ6XFxzKihbXixcXHNdKyl8dGltZW91dDpcXHMqKFxcZCspL2cpO1xyXG4gIGZvciAoY29uc3QgZiBvZiBwYWNrRnVuY3Rpb25zKSB7XHJcbiAgICBjb25zdCB0ZXN0cyA9IGYub3B0aW9uc1sndGVzdCddO1xyXG4gICAgY29uc3QgZGVtbyA9IGYub3B0aW9uc1snZGVtb1BhdGgnXTtcclxuICAgIGlmICgodGVzdHMgJiYgQXJyYXkuaXNBcnJheSh0ZXN0cykgJiYgdGVzdHMubGVuZ3RoKSkge1xyXG4gICAgICBmb3IgKGxldCBpID0gMDsgaSA8IHRlc3RzLmxlbmd0aDsgaSsrKSB7XHJcbiAgICAgICAgY29uc3QgcmVzID0gKHRlc3RzW2ldIGFzIHN0cmluZykubWF0Y2hBbGwocmVnKTtcclxuICAgICAgICBjb25zdCBtYXA6IHsgc2tpcD86IHN0cmluZywgd2FpdD86IG51bWJlciwgY2F0Pzogc3RyaW5nLCB0aW1lb3V0PzogbnVtYmVyLCBiZW5jaG1hcmtUaW1lb3V0PzogbnVtYmVyIH0gPSB7fTtcclxuICAgICAgICBBcnJheS5mcm9tKHJlcykuZm9yRWFjaCgoYXJyKSA9PiB7XHJcbiAgICAgICAgICBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ3NraXAnKSkgbWFwWydza2lwJ10gPSBhcnJbMV07XHJcbiAgICAgICAgICBlbHNlIGlmIChhcnJbMF0uc3RhcnRzV2l0aCgnd2FpdCcpKSBtYXBbJ3dhaXQnXSA9IHBhcnNlSW50KGFyclsyXSk7XHJcbiAgICAgICAgICBlbHNlIGlmIChhcnJbMF0uc3RhcnRzV2l0aCgnY2F0JykpIG1hcFsnY2F0J10gPSBhcnJbM107XHJcbiAgICAgICAgICBlbHNlIGlmIChhcnJbMF0uc3RhcnRzV2l0aCgndGltZW91dCcpKSBtYXBbJ3RpbWVvdXQnXSA9IHBhcnNlSW50KGFycls0XSk7XHJcbiAgICAgICAgfSk7XHJcbiAgICAgICAgY29uc3QgdGVzdCA9IG5ldyBUZXN0KG1hcC5jYXQgPz8gYXV0b1Rlc3RzQ2F0TmFtZSwgdGVzdHMubGVuZ3RoID09PSAxID8gZi5uYW1lIDogYCR7Zi5uYW1lfSAke2kgKyAxfWAsIGFzeW5jICgpID0+IHtcclxuICAgICAgICAgIGNvbnN0IHJlcyA9IGF3YWl0IGdyb2suZnVuY3Rpb25zLmV2YWwoYWRkTmFtZXNwYWNlKHRlc3RzW2ldLCBmKSk7XHJcbiAgICAgICAgICBpZiAobWFwLndhaXQpIGF3YWl0IGRlbGF5KG1hcC53YWl0KTtcclxuICAgICAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBuby10aHJvdy1saXRlcmFsXHJcbiAgICAgICAgICBpZiAodHlwZW9mIHJlcyA9PT0gJ2Jvb2xlYW4nICYmICFyZXMpIHRocm93IGBGYWlsZWQ6ICR7dGVzdHNbaV19LCBleHBlY3RlZCB0cnVlLCBnb3QgJHtyZXN9YDtcclxuICAgICAgICB9LCB7IHNraXBSZWFzb246IG1hcC5za2lwLCB0aW1lb3V0OiBERy5UZXN0LmlzSW5CZW5jaG1hcmsgPyBtYXAuYmVuY2htYXJrVGltZW91dCA/PyBCRU5DSE1BUktfVElNRU9VVCA6IG1hcC50aW1lb3V0ID8/IFNUQU5EQVJUX1RJTUVPVVQgfSk7XHJcbiAgICAgICAgaWYgKG1hcC5jYXQpIHtcclxuICAgICAgICAgIGNvbnN0IGNhdDogc3RyaW5nID0gbWFwLmNhdDtcclxuICAgICAgICAgIGlmIChtb2R1bGVUZXN0c1tjYXRdID09PSB1bmRlZmluZWQpXHJcbiAgICAgICAgICAgIG1vZHVsZVRlc3RzW2NhdF0gPSB7IHRlc3RzOiBbXSwgY2xlYXI6IHRydWUgfTtcclxuXHJcbiAgICAgICAgICAvLyBvbmx5IGJlZm9yZS9hZnRlciBjYW4gYmUgZGVmaW5lZCBpbiB0cyBmaWxlcyB0ZXN0cyB1bmRlciB0aGUgY2F0ZWdvcnlcclxuICAgICAgICAgIGlmICghbW9kdWxlVGVzdHNbY2F0XS50ZXN0cylcclxuICAgICAgICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cyA9IFtdO1xyXG4gICAgICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cy5wdXNoKHRlc3QpO1xyXG4gICAgICAgIH1cclxuICAgICAgICBlbHNlXHJcbiAgICAgICAgICBtb2R1bGVBdXRvVGVzdHMucHVzaCh0ZXN0KTtcclxuICAgICAgfVxyXG4gICAgfVxyXG4gICAgaWYgKGRlbW8pIHtcclxuICAgICAgY29uc3Qgd2FpdCA9IGYub3B0aW9uc1snZGVtb1dhaXQnXSA/IHBhcnNlSW50KGYub3B0aW9uc1snZGVtb1dhaXQnXSkgOiB1bmRlZmluZWQ7XHJcbiAgICAgIGNvbnN0IHRlc3QgPSBuZXcgVGVzdChkZW1vQ2F0TmFtZSwgZi5mcmllbmRseU5hbWUsIGFzeW5jICgpID0+IHtcclxuICAgICAgICBhd2FpdCBkZWxheSgzMDApO1xyXG4gICAgICAgIGdyb2suc2hlbGwuY2xlYXJMYXN0RXJyb3IoKTtcclxuICAgICAgICBhd2FpdCBmLmFwcGx5KCk7XHJcbiAgICAgICAgYXdhaXQgZGVsYXkod2FpdCA/IHdhaXQgOiAyMDAwKTtcclxuICAgICAgICBjb25zdCB1bmhhbmRsZWQgPSBhd2FpdCBncm9rLnNoZWxsLmxhc3RFcnJvcjtcclxuICAgICAgICBpZiAodW5oYW5kbGVkKVxyXG4gICAgICAgICAgdGhyb3cgbmV3IEVycm9yKHVuaGFuZGxlZCk7XHJcbiAgICAgIH0sIHsgc2tpcFJlYXNvbjogZi5vcHRpb25zWydkZW1vU2tpcCddIH0pO1xyXG4gICAgICBtb2R1bGVEZW1vLnB1c2godGVzdCk7XHJcbiAgICB9XHJcbiAgICBpZiAoZi5oYXNUYWcoJ3NlbVR5cGVEZXRlY3RvcicpKSB7XHJcbiAgICAgIGxldCBkZXRlY3RvcnNUZXN0RGF0YSA9IHRlc3REYXRhO1xyXG4gICAgICBpZiAoZi5vcHRpb25zWyd0ZXN0RGF0YSddKSB7XHJcbiAgICAgICAgZGV0ZWN0b3JzVGVzdERhdGEgPSBhd2FpdCBncm9rLmRhdGEuZmlsZXMub3BlblRhYmxlKGBTeXN0ZW06QXBwRGF0YS8ke3BhY2thZ2VfLm5xTmFtZX0vJHtmLm9wdGlvbnNbJ3Rlc3REYXRhJ119YCk7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGNvbnN0IHRlc3QgPSBuZXcgVGVzdChkZXRlY3RvcnNDYXROYW1lLCBmLmZyaWVuZGx5TmFtZSwgYXN5bmMgKCkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGFyciA9IFtdO1xyXG4gICAgICAgIGNvbnNvbGUubG9nKGBTeXN0ZW06QXBwRGF0YS8ke3BhY2thZ2VfLm5xTmFtZX0vJHtmLm9wdGlvbnNbJ3Rlc3REYXRhJ119YCk7XHJcblxyXG4gICAgICAgIGZvciAoY29uc3QgY29sIG9mIGRldGVjdG9yc1Rlc3REYXRhLmNsb25lKCkuY29sdW1ucykge1xyXG4gICAgICAgICAgY29uc3QgcmVzID0gYXdhaXQgZi5hcHBseShbY29sXSk7XHJcbiAgICAgICAgICBhcnIucHVzaChyZXMgfHwgY29sLnNlbVR5cGUpO1xyXG4gICAgICAgIH1cclxuICAgICAgICBjb25zdCByZXNBcnIgPSBhcnIuZmlsdGVyKChpKSA9PiBpKTtcclxuICAgICAgICBleHBlY3QocmVzQXJyLmxlbmd0aCwgMSk7XHJcblxyXG4gICAgICAgIGlmIChmLm9wdGlvbnNbJ3Rlc3REYXRhQ29sdW1uTmFtZSddKVxyXG4gICAgICAgICAgZXhwZWN0KHJlc0FyclswXSwgZi5vcHRpb25zWyd0ZXN0RGF0YUNvbHVtbk5hbWUnXSk7XHJcblxyXG4gICAgICB9LCB7IHNraXBSZWFzb246IGYub3B0aW9uc1snc2tpcFRlc3QnXSB9KTtcclxuICAgICAgbW9kdWxlRGV0ZWN0b3JzLnB1c2godGVzdCk7XHJcbiAgICB9XHJcbiAgfVxyXG4gIHdhc1JlZ2lzdGVyZWRbcGFja2FnZUlkXSA9IHRydWU7XHJcbiAgaWYgKG1vZHVsZUF1dG9UZXN0cy5sZW5ndGggPiAwKVxyXG4gICAgbW9kdWxlVGVzdHNbYXV0b1Rlc3RzQ2F0TmFtZV0gPSB7IHRlc3RzOiBtb2R1bGVBdXRvVGVzdHMsIGNsZWFyOiB0cnVlIH07XHJcbiAgaWYgKG1vZHVsZURlbW8ubGVuZ3RoID4gMClcclxuICAgIG1vZHVsZVRlc3RzW2RlbW9DYXROYW1lXSA9IHsgdGVzdHM6IG1vZHVsZURlbW8sIGNsZWFyOiB0cnVlIH07XHJcbiAgaWYgKG1vZHVsZURldGVjdG9ycy5sZW5ndGggPiAwKVxyXG4gICAgbW9kdWxlVGVzdHNbZGV0ZWN0b3JzQ2F0TmFtZV0gPSB7IHRlc3RzOiBtb2R1bGVEZXRlY3RvcnMsIGNsZWFyOiBmYWxzZSB9O1xyXG59XHJcblxyXG5mdW5jdGlvbiByZWRlZmluZUNvbnNvbGUoKTogYW55W10ge1xyXG4gIGNvbnN0IGxvZ3M6IGFueVtdID0gW107XHJcbiAgY29uc29sZS5sb2cgPSAoLi4uYXJncykgPT4ge1xyXG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xyXG4gICAgc3RkTG9nKC4uLmFyZ3MpO1xyXG4gIH07XHJcbiAgY29uc29sZS5pbmZvID0gKC4uLmFyZ3MpID0+IHtcclxuICAgIGxvZ3MucHVzaCguLi5hcmdzKTtcclxuICAgIHN0ZEluZm8oLi4uYXJncyk7XHJcbiAgfTtcclxuICBjb25zb2xlLndhcm4gPSAoLi4uYXJncykgPT4ge1xyXG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xyXG4gICAgc3RkV2FybiguLi5hcmdzKTtcclxuICB9O1xyXG4gIGNvbnNvbGUuZXJyb3IgPSAoLi4uYXJncykgPT4ge1xyXG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xyXG4gICAgc3RkRXJyb3IoLi4uYXJncyk7XHJcbiAgfTtcclxuICByZXR1cm4gbG9ncztcclxufVxyXG5cclxuZnVuY3Rpb24gcmVzZXRDb25zb2xlKCk6IHZvaWQge1xyXG4gIGNvbnNvbGUubG9nID0gc3RkTG9nO1xyXG4gIGNvbnNvbGUuaW5mbyA9IHN0ZEluZm87XHJcbiAgY29uc29sZS53YXJuID0gc3RkV2FybjtcclxuICBjb25zb2xlLmVycm9yID0gc3RkRXJyb3I7XHJcbn1cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBydW5UZXN0cyhvcHRpb25zPzogVGVzdEV4ZWN1dGlvbk9wdGlvbnMpIDogUHJvbWlzZTxUZXN0UmVzdWx0RXh0ZW5kZWRbXT57XHJcblxyXG4gIGNvbnN0IHBhY2thZ2VfOiBfREcuUGFja2FnZSA9IG9wdGlvbnM/Lm5vZGVPcHRpb25zID8gb3B0aW9ucy5ub2RlT3B0aW9ucy5wYWNrYWdlIDogZ3Jvay5mdW5jdGlvbnMuZ2V0Q3VycmVudENhbGwoKS5mdW5jLnBhY2thZ2U7XHJcbiAgaWYgKCFwYWNrYWdlXylcclxuICAgIHRocm93IG5ldyBFcnJvcignQ2FuXFwndCBydW4gdGVzdHMgb3V0c2lkZSBvZiB0aGUgcGFja2FnZScpO1xyXG4gIGNvbnN0IG1hdGNoID0gcGFja2FnZV8ucGFja2FnZU93bmVyPy5tYXRjaCgvPChbXj5dKik+Lyk7XHJcbiAgY29uc3QgcGFja2FnZU93bmVyID0gbWF0Y2ggPyBtYXRjaFsxXSA6ICcnO1xyXG4gIGlmIChwYWNrYWdlXyAhPSB1bmRlZmluZWQpXHJcbiAgICBhd2FpdCBpbml0QXV0b1Rlc3RzKHBhY2thZ2VfKTtcclxuICBjb25zdCByZXN1bHRzOlRlc3RSZXN1bHRFeHRlbmRlZFtdID0gW107XHJcbiAgY29uc29sZS5sb2coYFJ1bm5pbmcgdGVzdHNgKTtcclxuICBjb25zb2xlLmxvZyhvcHRpb25zKTtcclxuICBvcHRpb25zID8/PSB7fTtcclxuICBvcHRpb25zIS50ZXN0Q29udGV4dCA/Pz0gbmV3IFRlc3RDb250ZXh0KCk7XHJcbiAgZ3Jvay5zaGVsbC5jbGVhckxhc3RFcnJvcigpO1xyXG4gIGNvbnN0IGxvZ3MgPSByZWRlZmluZUNvbnNvbGUoKTtcclxuXHJcbiAgYXdhaXQgaW52b2tlVGVzdHModGVzdHMsIG9wdGlvbnMpO1xyXG5cclxuICBmb3IgKGxldCByIG9mIHJlc3VsdHMpIHtcclxuICAgIHIucmVzdWx0ID0gci5yZXN1bHQudG9TdHJpbmcoKS5yZXBsYWNlKC9cIi9nLCAnXFwnJyk7XHJcbiAgICBpZiAoci5sb2dzICE9IHVuZGVmaW5lZClcclxuICAgICAgci5sb2dzID0gci5sb2dzIS50b1N0cmluZygpLnJlcGxhY2UoL1wiL2csICdcXCcnKTtcclxuICB9XHJcbiAgcmV0dXJuIHJlc3VsdHM7XHJcblxyXG4gIGFzeW5jIGZ1bmN0aW9uIGludm9rZUNhdGVnb3J5TWV0aG9kKG1ldGhvZDogKCgpID0+IFByb21pc2U8dm9pZD4pIHwgdW5kZWZpbmVkLCBjYXRlZ29yeTogc3RyaW5nKTogUHJvbWlzZTxzdHJpbmcgfCB1bmRlZmluZWQ+IHtcclxuICAgIGxldCBpbnZva2F0aW9uUmVzdWx0ID0gdW5kZWZpbmVkO1xyXG4gICAgdHJ5IHtcclxuICAgICAgaWYgKG1ldGhvZCAhPT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgYXdhaXQgdGltZW91dChhc3luYyAoKSA9PiB7XHJcbiAgICAgICAgICBhd2FpdCBtZXRob2QoKTtcclxuICAgICAgICB9LCAxMDAwMDAsIGBiZWZvcmUgJHtjYXRlZ29yeX06IHRpbWVvdXQgZXJyb3JgKTtcclxuICAgICAgfVxyXG4gICAgfSBjYXRjaCAoeDogYW55KSB7XHJcbiAgICAgIGludm9rYXRpb25SZXN1bHQgPSBhd2FpdCBnZXRSZXN1bHQoeCk7XHJcbiAgICB9XHJcbiAgICByZXR1cm4gaW52b2thdGlvblJlc3VsdFxyXG4gIH1cclxuXHJcbiAgYXN5bmMgZnVuY3Rpb24gaW52b2tlVGVzdHNJbkNhdGVnb3J5KGNhdGVnb3J5OiBDYXRlZ29yeSwgb3B0aW9uczogVGVzdEV4ZWN1dGlvbk9wdGlvbnMsIGlzVGFyZ2V0Q2F0ZWdvcnk6IGJvb2xlYW4pOiBQcm9taXNlPFRlc3RSZXN1bHRFeHRlbmRlZFtdPiB7XHJcbiAgICBsZXQgdCA9IGNhdGVnb3J5LnRlc3RzID8/IFtdO1xyXG4gICAgY29uc3QgcmVzIDogVGVzdFJlc3VsdEV4dGVuZGVkW10gPSBbXTtcclxuICAgIC8vIGxldCBtZW1vcnlVc2FnZUJlZm9yZSA9ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemU7XHJcbiAgICBjb25zdCB3aWRnZXRzQmVmb3JlID0gZ2V0V2lkZ2V0c0NvdW50U2FmZSgpO1xyXG5cclxuICAgIGlmIChjYXRlZ29yeS5jbGVhcikge1xyXG4gICAgICAgIGxldCBza2lwcGluZ1Rlc3RzID0gaXNUYXJnZXRDYXRlZ29yeSAmJiBvcHRpb25zLnNraXBUb1Rlc3QgIT0gdW5kZWZpbmVkO1xyXG4gICAgICBmb3IgKGxldCBpID0gMDsgaSA8IHQubGVuZ3RoOyBpKyspIHtcclxuXHJcbiAgICAgICAgaWYgKHRbaV0ub3B0aW9ucykge1xyXG4gICAgICAgICAgaWYgKHRbaV0ub3B0aW9ucz8uYmVuY2htYXJrID09PSB1bmRlZmluZWQpIHtcclxuICAgICAgICAgICAgaWYgKCF0W2ldLm9wdGlvbnMpXHJcbiAgICAgICAgICAgICAgdFtpXS5vcHRpb25zID0ge31cclxuICAgICAgICAgICAgdFtpXS5vcHRpb25zIS5iZW5jaG1hcmsgPSBjYXRlZ29yeS5iZW5jaG1hcmtzID8/IGZhbHNlO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICAgICBsZXQgdGVzdCA9IHRbaV07XHJcbiAgICAgICAgaWYgKG9wdGlvbnMudGVzdClcclxuICAgICAgICAgIGlmIChvcHRpb25zLnRlc3QudG9Mb3dlckNhc2UoKSAhPT0gdGVzdC5uYW1lLnRvTG93ZXJDYXNlKCkpXHJcbiAgICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgIGlmIChza2lwcGluZ1Rlc3RzKSB7XHJcbiAgICAgICAgICBpZiAob3B0aW9ucz8uc2tpcFRvVGVzdCAhPSB1bmRlZmluZWQgJiYgdGVzdC5uYW1lLnRvTG93ZXJDYXNlKCkudHJpbSgpID09PSBvcHRpb25zPy5za2lwVG9UZXN0LnRvTG93ZXJDYXNlKCkudHJpbSgpKSB7XHJcbiAgICAgICAgICAgIC8vIEZvdW5kIHRoZSB0YXJnZXQgdGVzdCwgc3RvcCBza2lwcGluZyBhZnRlciB0aGlzIG9uZVxyXG4gICAgICAgICAgICBza2lwcGluZ1Rlc3RzID0gZmFsc2U7XHJcbiAgICAgICAgICB9IGVsc2VcclxuICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgIH1cclxuICAgICAgICBpZiAodGVzdD8ub3B0aW9ucykge1xyXG4gICAgICAgICAgdGVzdC5vcHRpb25zLm93bmVyID0gdFtpXS5vcHRpb25zPy5vd25lciA/PyBjYXRlZ29yeT8ub3duZXIgPz8gcGFja2FnZU93bmVyID8/ICcnO1xyXG4gICAgICAgIH1cclxuICAgICAgICAvLyBsZXQgaXNHQkVuYWJsZSA9ICh3aW5kb3cgYXMgYW55KS5nYyAmJiB0ZXN0Lm9wdGlvbnM/LnNraXBSZWFzb24gPT0gdW5kZWZpbmVkO1xyXG4gICAgICAgIC8vIGNvbnNvbGUubG9nKGAqKioqKioqKiR7aXNHQkVuYWJsZX1gKTtcclxuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcclxuICAgICAgICAvLyAgIGF3YWl0ICh3aW5kb3cgYXMgYW55KS5nYygpO1xyXG4gICAgICAgIC8vIG1lbW9yeVVzYWdlQmVmb3JlID0gKHdpbmRvdz8ucGVyZm9ybWFuY2UgYXMgYW55KT8ubWVtb3J5Py51c2VkSlNIZWFwU2l6ZTtcclxuICAgICAgICBsZXQgdGVzdFJ1biA9IGF3YWl0IGV4ZWNUZXN0KFxyXG4gICAgICAgICAgICB0ZXN0LFxyXG4gICAgICAgICAgICBvcHRpb25zPy50ZXN0LFxyXG4gICAgICAgICAgICBsb2dzLCBERy5UZXN0LmlzSW5CZW5jaG1hcmsgPyB0W2ldLm9wdGlvbnM/LmJlbmNobWFya1RpbWVvdXQgPz8gQkVOQ0hNQVJLX1RJTUVPVVQgOiB0W2ldLm9wdGlvbnM/LnRpbWVvdXQgPz8gU1RBTkRBUlRfVElNRU9VVCxcclxuICAgICAgICAgICAgcGFja2FnZV8ubmFtZSxcclxuICAgICAgICAgICAgb3B0aW9ucy52ZXJib3NlXHJcbiAgICAgICAgKTtcclxuXHJcbiAgICAgICAgLy8gaWYgKGlzR0JFbmFibGUpXHJcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcclxuICAgICAgICBpZiAodGVzdFJ1bikge1xyXG4gICAgICAgICAgcmVzLnB1c2goeyAuLi50ZXN0UnVuLCAgd2lkZ2V0c0RpZmZlcmVuY2U6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XHJcbiAgICAgICAgICAvLyBSZXR1cm4gZWFybHkgaWYgcmV0dXJuT25GYWlsIGlzIHNldCBhbmQgdGVzdCBmYWlsZWQgKGJ1dCBpZ25vcmUgZmFpbHVyZSBmb3IgdGhlIHNraXBUb1Rlc3QgdGVzdCBpdHNlbGYpXHJcbiAgICAgICAgICBpZiAob3B0aW9ucy5yZXR1cm5PbkZhaWwgJiYgb3B0aW9ucy5za2lwVG9UZXN0ICE9PSB0ZXN0Lm5hbWUgJiYgIXRlc3RSdW4uc3VjY2VzcyAmJiAhdGVzdFJ1bi5za2lwcGVkKVxyXG4gICAgICAgICAgICByZXR1cm4gcmVzO1xyXG4gICAgICAgIH1cclxuICAgICAgICAvLyByZXMucHVzaCh7IC4uLnRlc3RSdW4sIG1lbW9yeURlbHRhOiAod2luZG93Py5wZXJmb3JtYW5jZSBhcyBhbnkpPy5tZW1vcnk/LnVzZWRKU0hlYXBTaXplIC0gbWVtb3J5VXNhZ2VCZWZvcmUsIHdpZGdldHNEZWx0YTogZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIC0gd2lkZ2V0c0JlZm9yZSB9KTtcclxuXHJcbiAgICAgICAgaWYgKCFvcHRpb25zLm5vZGVPcHRpb25zKSB7XHJcbiAgICAgICAgICBncm9rLnNoZWxsLmNsb3NlQWxsKCk7XHJcbiAgICAgICAgICBERy5CYWxsb29uLmNsb3NlQWxsKCk7XHJcbiAgICAgICAgfVxyXG4gICAgICB9XHJcbiAgICB9IGVsc2Uge1xyXG4gICAgICBsZXQgc2tpcHBpbmdUZXN0cyA9IGlzVGFyZ2V0Q2F0ZWdvcnkgJiYgb3B0aW9ucy5za2lwVG9UZXN0ICE9IHVuZGVmaW5lZDtcclxuICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCB0Lmxlbmd0aDsgaSsrKSB7XHJcbiAgICAgICAgbGV0IHRlc3QgPSB0W2ldO1xyXG4gICAgICAgIGlmIChvcHRpb25zLnRlc3QpXHJcbiAgICAgICAgICBpZiAob3B0aW9ucy50ZXN0LnRvTG93ZXJDYXNlKCkgIT09IHRlc3QubmFtZS50b0xvd2VyQ2FzZSgpKVxyXG4gICAgICAgICAgICBjb250aW51ZTtcclxuICAgICAgICBpZiAoc2tpcHBpbmdUZXN0cykge1xyXG4gICAgICAgICAgaWYgKG9wdGlvbnM/LnNraXBUb1Rlc3QgIT0gdW5kZWZpbmVkICYmIHRlc3QubmFtZS50b0xvd2VyQ2FzZSgpLnRyaW0oKSA9PT0gb3B0aW9ucz8uc2tpcFRvVGVzdC50b0xvd2VyQ2FzZSgpLnRyaW0oKSkge1xyXG4gICAgICAgICAgICAvLyBGb3VuZCB0aGUgdGFyZ2V0IHRlc3QsIHN0b3Agc2tpcHBpbmcgYWZ0ZXIgdGhpcyBvbmVcclxuICAgICAgICAgICAgc2tpcHBpbmdUZXN0cyA9IGZhbHNlO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgY29udGludWU7ICAvLyBTa2lwIHRoaXMgdGVzdCAoaW5jbHVkaW5nIHRoZSB0YXJnZXQpXHJcbiAgICAgICAgfVxyXG5cclxuICAgICAgICBpZiAodGVzdD8ub3B0aW9ucykge1xyXG4gICAgICAgICAgdGVzdC5vcHRpb25zLm93bmVyID0gdFtpXS5vcHRpb25zPy5vd25lciA/PyBjYXRlZ29yeT8ub3duZXIgPz8gcGFja2FnZU93bmVyID8/ICcnO1xyXG4gICAgICAgIH1cclxuICAgICAgICAvLyBsZXQgaXNHQkVuYWJsZSA9ICh3aW5kb3cgYXMgYW55KS5nYyAmJiB0ZXN0Lm9wdGlvbnM/LnNraXBSZWFzb24gPT0gdW5kZWZpbmVkO1xyXG4gICAgICAgIC8vIGNvbnNvbGUubG9nKGAqKioqKioqKiR7aXNHQkVuYWJsZX1gKTtcclxuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcclxuICAgICAgICAvLyAgIGF3YWl0ICh3aW5kb3cgYXMgYW55KS5nYygpO1xyXG4gICAgICAgIC8vIG1lbW9yeVVzYWdlQmVmb3JlID0gKHdpbmRvdz8ucGVyZm9ybWFuY2UgYXMgYW55KT8ubWVtb3J5Py51c2VkSlNIZWFwU2l6ZTtcclxuICAgICAgICBsZXQgdGVzdFJ1biA9IGF3YWl0IGV4ZWNUZXN0KFxyXG4gICAgICAgICAgICB0ZXN0LFxyXG4gICAgICAgICAgICBvcHRpb25zPy50ZXN0LFxyXG4gICAgICAgICAgICBsb2dzLFxyXG4gICAgICAgICAgICBERy5UZXN0LmlzSW5CZW5jaG1hcmsgPyB0W2ldLm9wdGlvbnM/LmJlbmNobWFya1RpbWVvdXQgPz8gQkVOQ0hNQVJLX1RJTUVPVVQgOiB0W2ldLm9wdGlvbnM/LnRpbWVvdXQsXHJcbiAgICAgICAgICAgIHBhY2thZ2VfLm5hbWUsXHJcbiAgICAgICAgICAgIG9wdGlvbnMudmVyYm9zZVxyXG4gICAgICAgICk7XHJcblxyXG4gICAgICAgIC8vIGlmIChpc0dCRW5hYmxlKVxyXG4gICAgICAgIC8vICAgYXdhaXQgKHdpbmRvdyBhcyBhbnkpLmdjKCk7XHJcblxyXG4gICAgICAgIGlmICh0ZXN0UnVuKSB7XHJcbiAgICAgICAgICByZXMucHVzaCh7IC4uLnRlc3RSdW4sIHdpZGdldHNEaWZmZXJlbmNlOiBnZXRXaWRnZXRzQ291bnRTYWZlKCkgLSB3aWRnZXRzQmVmb3JlIH0pO1xyXG4gICAgICAgICAgLy8gUmV0dXJuIGVhcmx5IGlmIHJldHVybk9uRmFpbCBpcyBzZXQgYW5kIHRlc3QgZmFpbGVkIChidXQgaWdub3JlIGZhaWx1cmUgZm9yIHRoZSBza2lwVG9UZXN0IHRlc3QgaXRzZWxmKVxyXG4gICAgICAgICAgaWYgKG9wdGlvbnMucmV0dXJuT25GYWlsICYmIG9wdGlvbnMuc2tpcFRvVGVzdCAhPT0gdGVzdC5uYW1lICYmICF0ZXN0UnVuLnN1Y2Nlc3MgJiYgIXRlc3RSdW4uc2tpcHBlZClcclxuICAgICAgICAgICAgcmV0dXJuIHJlcztcclxuICAgICAgICB9XHJcbiAgICAgICAgLy8gcmVzLnB1c2goeyAuLi50ZXN0UnVuLCBtZW1vcnlEZWx0YTogKHdpbmRvdz8ucGVyZm9ybWFuY2UgYXMgYW55KT8ubWVtb3J5Py51c2VkSlNIZWFwU2l6ZSAtIG1lbW9yeVVzYWdlQmVmb3JlLCB3aWRnZXRzRGlmZmVyZW5jZTogZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIC0gd2lkZ2V0c0JlZm9yZSB9KTtcclxuXHJcbiAgICAgIH1cclxuICAgIH1cclxuICAgIHJldHVybiByZXM7XHJcbiAgfVxyXG5cclxuICBmdW5jdGlvbiBnZXRXaWRnZXRzQ291bnRTYWZlKCkge1xyXG4gICAgaWYgKHR5cGVvZiBwcm9jZXNzICE9PSAndW5kZWZpbmVkJylcclxuICAgICAgcmV0dXJuIDA7XHJcbiAgICBsZXQgbGVuZ3RoID0gLTE7XHJcbiAgICB0cnkge1xyXG4gICAgICBsZW5ndGggPSBERy5XaWRnZXQuZ2V0QWxsKCkubGVuZ3RoO1xyXG4gICAgfSBjYXRjaCAoZTogYW55KSB7XHJcbiAgICAgIGNvbnNvbGUud2FybihlLm1lc3NhZ2UgPz8gZSk7XHJcbiAgICB9XHJcbiAgICByZXR1cm4gbGVuZ3RoO1xyXG4gIH1cclxuXHJcbiAgYXN5bmMgZnVuY3Rpb24gaW52b2tlVGVzdHMoY2F0ZWdvcmllc1RvSW52b2tlOiB7IFtrZXk6IHN0cmluZ106IENhdGVnb3J5IH0sIG9wdGlvbnM6IFRlc3RFeGVjdXRpb25PcHRpb25zKSB7XHJcbiAgICB0cnkge1xyXG4gICAgICBsZXQgc2tpcHBpbmdDYXRlZ29yaWVzID0gb3B0aW9ucz8uc2tpcFRvQ2F0ZWdvcnkgIT0gdW5kZWZpbmVkO1xyXG4gICAgICBsZXQgaXNUYXJnZXRDYXRlZ29yeSA9IGZhbHNlO1xyXG4gICAgICBmb3IgKGNvbnN0IFtrZXksIHZhbHVlXSBvZiBPYmplY3QuZW50cmllcyhjYXRlZ29yaWVzVG9JbnZva2UpKSB7XHJcbiAgICAgICAgICBpZiAob3B0aW9ucy5leGNsdWRlPy5zb21lKChjKSA9PiBrZXkuc3RhcnRzV2l0aChjKSkpXHJcbiAgICAgICAgICAgICAgY29udGludWU7XHJcbiAgICAgICAgICBpZiAob3B0aW9ucz8uY2F0ZWdvcnkgIT0gbnVsbCAmJiAha2V5LnRvTG93ZXJDYXNlKCkuc3RhcnRzV2l0aChgJHtvcHRpb25zPy5jYXRlZ29yeS50b0xvd2VyQ2FzZSgpLnRyaW0oKX0gOmApICYmXHJcbiAgICAgICAgICAgICAga2V5LnRvTG93ZXJDYXNlKCkudHJpbSgpICE9PSBvcHRpb25zPy5jYXRlZ29yeS50b0xvd2VyQ2FzZSgpLnRyaW0oKSlcclxuICAgICAgICAgICAgICBjb250aW51ZTtcclxuXHJcbiAgICAgICAgICBpZiAoc2tpcHBpbmdDYXRlZ29yaWVzKSB7XHJcbiAgICAgICAgICAgICAgaWYgKGlzVGFyZ2V0Q2F0ZWdvcnkpXHJcbiAgICAgICAgICAgICAgICAgIHNraXBwaW5nQ2F0ZWdvcmllcyA9IGZhbHNlO1xyXG4gICAgICAgICAgICAgIGVsc2Uge1xyXG4gICAgICAgICAgICAgICAgICBpZiAob3B0aW9ucz8uc2tpcFRvQ2F0ZWdvcnkgIT0gbnVsbCAmJiBrZXkudG9Mb3dlckNhc2UoKS50cmltKCkgPT09IG9wdGlvbnM/LnNraXBUb0NhdGVnb3J5LnRvTG93ZXJDYXNlKCkudHJpbSgpKSB7XHJcbiAgICAgICAgICAgICAgICAgICAgICBpc1RhcmdldENhdGVnb3J5ID0gdHJ1ZTtcclxuICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcclxuICAgICAgICAgICAgICAgICAgICAgIC8vIEhhdmVuJ3QgZm91bmQgdGhlIHRhcmdldCBjYXRlZ29yeSB5ZXQsIGtlZXAgc2tpcHBpbmdcclxuICAgICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgICAgICAgICAgICB9XHJcbiAgICAgICAgICAgICAgfVxyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFN0YXJ0ZWQge3ske2tleX19fWApO1xyXG4gICAgICAgICAgLy9AdHMtaWdub3JlXHJcbiAgICAgICAgICBjb25zdCBza2lwcGVkID0gdmFsdWUudGVzdHM/LmV2ZXJ5KCh0OiBUZXN0KSA9PiB0Lm9wdGlvbnM/LnNraXBSZWFzb24pO1xyXG4gICAgICAgICAgaWYgKCFza2lwcGVkKVxyXG4gICAgICAgICAgICAgIHZhbHVlLmJlZm9yZVN0YXR1cyA9IGF3YWl0IGludm9rZUNhdGVnb3J5TWV0aG9kKHZhbHVlLmJlZm9yZSwga2V5KTtcclxuXHJcbiAgICAgICAgICBsZXQgdCA9IHZhbHVlLnRlc3RzID8/IFtdO1xyXG5cclxuICAgICAgICAgIGlmIChvcHRpb25zLnN0cmVzc1Rlc3QpIHtcclxuICAgICAgICAgICAgICB0ID0gdC5maWx0ZXIoKGUpID0+IGUub3B0aW9ucz8uc3RyZXNzVGVzdCk7XHJcbiAgICAgICAgICAgICAgdCA9IHNodWZmbGUodCk7XHJcbiAgICAgICAgICB9XHJcblxyXG4gICAgICAgICAgaWYgKChvcHRpb25zLnRhZ3M/Lmxlbmd0aCA/PyAwKSA+IDApIHtcclxuICAgICAgICAgICAgICB0ID0gdC5maWx0ZXIoKGUpID0+XHJcbiAgICAgICAgICAgICAgICAgIGUub3B0aW9ucz8udGFncz8uc29tZSh0YWcgPT4gKG9wdGlvbnM/LnRhZ3MgPz8gW10pLmluY2x1ZGVzKHRhZykpXHJcbiAgICAgICAgICAgICAgKTtcclxuICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICBsZXQgcmVzOiBUZXN0UmVzdWx0RXh0ZW5kZWRbXTtcclxuICAgICAgICAgIGlmICh2YWx1ZS5iZWZvcmVTdGF0dXMpIHtcclxuICAgICAgICAgICAgICByZXMgPSBBcnJheS5mcm9tKHQubWFwKCh0ZXN0RWxlbSkgPT4ge1xyXG4gICAgICAgICAgICAgICAgICByZXR1cm4ge1xyXG4gICAgICAgICAgICAgICAgICAgICAgZGF0ZTogbmV3IERhdGUoKS50b0lTT1N0cmluZygpLFxyXG4gICAgICAgICAgICAgICAgICAgICAgY2F0ZWdvcnk6IGtleSxcclxuICAgICAgICAgICAgICAgICAgICAgIG5hbWU6IHRlc3RFbGVtLm5hbWUsXHJcbiAgICAgICAgICAgICAgICAgICAgICBzdWNjZXNzOiBmYWxzZSxcclxuICAgICAgICAgICAgICAgICAgICAgIHJlc3VsdDogJ2JlZm9yZSgpIGZhaWxlZCcsXHJcbiAgICAgICAgICAgICAgICAgICAgICBtczogMCxcclxuICAgICAgICAgICAgICAgICAgICAgIHNraXBwZWQ6IGZhbHNlLFxyXG4gICAgICAgICAgICAgICAgICAgICAgbG9nczogJycsXHJcbiAgICAgICAgICAgICAgICAgICAgICBvd25lcjogcGFja2FnZU93bmVyLFxyXG4gICAgICAgICAgICAgICAgICAgICAgcGFja2FnZTogcGFja2FnZV8ubmFtZSxcclxuICAgICAgICAgICAgICAgICAgICAgIHdpZGdldHNEaWZmZXJlbmNlOiAwLFxyXG4gICAgICAgICAgICAgICAgICAgICAgZmxha2luZzogREcuVGVzdC5pc1JlcHJvZHVjaW5nXHJcbiAgICAgICAgICAgICAgICAgIH07XHJcbiAgICAgICAgICAgICAgfSkpO1xyXG4gICAgICAgICAgICAgIHJlcy5mb3JFYWNoKGFzeW5jICh0ZXN0KSA9PiBhd2FpdCBncm9rLnNoZWxsLnJlcG9ydFRlc3QoJ3BhY2thZ2UnLCB0ZXN0KSk7XHJcbiAgICAgICAgICB9IGVsc2VcclxuICAgICAgICAgICAgICByZXMgPSBhd2FpdCBpbnZva2VUZXN0c0luQ2F0ZWdvcnkodmFsdWUsIG9wdGlvbnMsIHNraXBwaW5nQ2F0ZWdvcmllcyk7XHJcbiAgICAgICAgICBjb25zdCBkYXRhOiBUZXN0UmVzdWx0RXh0ZW5kZWRbXSA9IHJlcy5maWx0ZXIoKGQpID0+IGQucmVzdWx0ICE9ICdza2lwcGVkJyk7XHJcblxyXG4gICAgICAgICAgaWYgKCFza2lwcGVkKVxyXG4gICAgICAgICAgICAgIHZhbHVlLmFmdGVyU3RhdHVzID0gYXdhaXQgaW52b2tlQ2F0ZWdvcnlNZXRob2QodmFsdWUuYWZ0ZXIsIGtleSk7XHJcblxyXG4gICAgICAgICAgLy8gQ2xlYXIgYWZ0ZXIgY2F0ZWdvcnlcclxuICAgICAgICAgIC8vIGdyb2suc2hlbGwuY2xvc2VBbGwoKTtcclxuICAgICAgICAgIC8vIERHLkJhbGxvb24uY2xvc2VBbGwoKTtcclxuICAgICAgICAgIGlmICh2YWx1ZS5hZnRlclN0YXR1cykge1xyXG4gICAgICAgICAgICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBDYXRlZ29yeSBhZnRlcigpIHt7JHtrZXl9fX0gZmFpbGVkYCk7XHJcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFJlc3VsdCBmb3Ige3ske2tleX19fSBhZnRlcjogJHt2YWx1ZS5hZnRlclN0YXR1c31gKTtcclxuICAgICAgICAgICAgICBkYXRhLnB1c2goe1xyXG4gICAgICAgICAgICAgICAgICBkYXRlOiBuZXcgRGF0ZSgpLnRvSVNPU3RyaW5nKCksXHJcbiAgICAgICAgICAgICAgICAgIGNhdGVnb3J5OiBrZXksXHJcbiAgICAgICAgICAgICAgICAgIG5hbWU6ICdhZnRlcicsXHJcbiAgICAgICAgICAgICAgICAgIHN1Y2Nlc3M6IGZhbHNlLFxyXG4gICAgICAgICAgICAgICAgICByZXN1bHQ6IHZhbHVlLmFmdGVyU3RhdHVzLFxyXG4gICAgICAgICAgICAgICAgICBtczogMCxcclxuICAgICAgICAgICAgICAgICAgc2tpcHBlZDogZmFsc2UsXHJcbiAgICAgICAgICAgICAgICAgIGxvZ3M6ICcnLFxyXG4gICAgICAgICAgICAgICAgICBvd25lcjogcGFja2FnZU93bmVyLFxyXG4gICAgICAgICAgICAgICAgICBwYWNrYWdlOiBwYWNrYWdlXy5uYW1lLFxyXG4gICAgICAgICAgICAgICAgICB3aWRnZXRzRGlmZmVyZW5jZTogMCxcclxuICAgICAgICAgICAgICAgICAgZmxha2luZzogREcuVGVzdC5pc1JlcHJvZHVjaW5nXHJcbiAgICAgICAgICAgICAgfSk7XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgICBpZiAodmFsdWUuYmVmb3JlU3RhdHVzKSB7XHJcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IENhdGVnb3J5IGJlZm9yZSgpIHt7JHtrZXl9fX0gZmFpbGVkYCk7XHJcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFJlc3VsdCBmb3Ige3ske2tleX19fSBiZWZvcmU6ICR7dmFsdWUuYmVmb3JlU3RhdHVzfWApO1xyXG4gICAgICAgICAgICAgIGRhdGEucHVzaCh7XHJcbiAgICAgICAgICAgICAgICAgIGRhdGU6IG5ldyBEYXRlKCkudG9JU09TdHJpbmcoKSxcclxuICAgICAgICAgICAgICAgICAgY2F0ZWdvcnk6IGtleSxcclxuICAgICAgICAgICAgICAgICAgbmFtZTogJ2JlZm9yZScsXHJcbiAgICAgICAgICAgICAgICAgIHN1Y2Nlc3M6IGZhbHNlLFxyXG4gICAgICAgICAgICAgICAgICByZXN1bHQ6IHZhbHVlLmJlZm9yZVN0YXR1cyxcclxuICAgICAgICAgICAgICAgICAgbXM6IDAsXHJcbiAgICAgICAgICAgICAgICAgIHNraXBwZWQ6IGZhbHNlLFxyXG4gICAgICAgICAgICAgICAgICBsb2dzOiAnJyxcclxuICAgICAgICAgICAgICAgICAgb3duZXI6IHBhY2thZ2VPd25lcixcclxuICAgICAgICAgICAgICAgICAgcGFja2FnZTogcGFja2FnZV8ubmFtZSxcclxuICAgICAgICAgICAgICAgICAgd2lkZ2V0c0RpZmZlcmVuY2U6IDAsXHJcbiAgICAgICAgICAgICAgICAgIGZsYWtpbmc6IERHLlRlc3QuaXNSZXByb2R1Y2luZ1xyXG4gICAgICAgICAgICAgIH0pO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgcmVzdWx0cy5wdXNoKC4uLmRhdGEpO1xyXG5cclxuICAgICAgICAgIC8vIElmIHJldHVybk9uRmFpbCBpcyBzZXQgYW5kIGEgdGVzdCBmYWlsZWQgKG90aGVyIHRoYW4gc2tpcFRvVGVzdCksIHN0b3AgcHJvY2Vzc2luZyBtb3JlIGNhdGVnb3JpZXNcclxuICAgICAgICAgIGlmIChvcHRpb25zLnJldHVybk9uRmFpbCAmJiBkYXRhLnNvbWUoKGQpID0+ICFkLnN1Y2Nlc3MgJiYgIWQuc2tpcHBlZCAmJiBkLm5hbWUgIT09IG9wdGlvbnMuc2tpcFRvVGVzdCkpXHJcbiAgICAgICAgICAgICAgYnJlYWs7XHJcbiAgICAgIH1cclxuICAgIH0gZmluYWxseSB7XHJcbiAgICAgIHJlc2V0Q29uc29sZSgpO1xyXG4gICAgfVxyXG4gICAgaWYgKG9wdGlvbnMudGVzdENvbnRleHQhLmNhdGNoVW5oYW5kbGVkICYmICghREcuVGVzdC5pc0luQmVuY2htYXJrKSkge1xyXG4gICAgICBhd2FpdCBkZWxheSgxMDAwKTtcclxuICAgICAgY29uc3QgZXJyb3IgPSBhd2FpdCBncm9rLnNoZWxsLmxhc3RFcnJvcjtcclxuICAgICAgaWYgKGVycm9yICE9IHVuZGVmaW5lZCkge1xyXG4gICAgICAgICAgY29uc3QgcGFyYW1zOiBhbnkgPSB7XHJcbiAgICAgICAgICAgICAgbG9nczogJycsXHJcbiAgICAgICAgICAgICAgZGF0ZTogbmV3IERhdGUoKS50b0lTT1N0cmluZygpLFxyXG4gICAgICAgICAgICAgIGNhdGVnb3J5OiAnVW5oYW5kbGVkIGV4Y2VwdGlvbnMnLFxyXG4gICAgICAgICAgICAgIG5hbWU6ICdFeGNlcHRpb24nLFxyXG4gICAgICAgICAgICAgIHJlc3VsdDogZXJyb3IgPz8gJycsXHJcbiAgICAgICAgICAgICAgc3VjY2VzczogIWVycm9yLFxyXG4gICAgICAgICAgICAgIG1zOiAwLFxyXG4gICAgICAgICAgICAgIHNraXBwZWQ6IGZhbHNlLFxyXG4gICAgICAgICAgICAgIG93bmVyOiBwYWNrYWdlT3duZXIgPz8gJycsXHJcbiAgICAgICAgICAgICAgJ3BhY2thZ2UnOiBwYWNrYWdlXy5uYW1lLFxyXG4gICAgICAgICAgICAgIHdpZGdldHNEaWZmZXJlbmNlOiAwXHJcbiAgICAgICAgICB9O1xyXG4gICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFVuaGFuZGxlZCBFeGNlcHRpb246ICR7ZXJyb3J9YCk7XHJcblxyXG4gICAgICAgICAgcmVzdWx0cy5wdXNoKHsuLi5wYXJhbXMsICdmbGFraW5nJzogREcuVGVzdC5pc1JlcHJvZHVjaW5nICYmICFlcnJvcn0pO1xyXG4gICAgICAgICAgKDxhbnk+cGFyYW1zKS5wYWNrYWdlID0gcGFja2FnZV8ubmFtZTtcclxuICAgICAgICAgIGF3YWl0IGdyb2suc2hlbGwucmVwb3J0VGVzdCgncGFja2FnZScsIHBhcmFtcyk7XHJcbiAgICAgIH1cclxuICAgIH1cclxuICB9XHJcbn1cclxuXHJcbmFzeW5jIGZ1bmN0aW9uIGdldFJlc3VsdCh4OiBhbnkpOiBQcm9taXNlPHN0cmluZz4ge1xyXG4gIHJldHVybiBgJHt4LnRvU3RyaW5nKCl9XFxuJHt4LnN0YWNrID8gKGF3YWl0IERHLkxvZ2dlci50cmFuc2xhdGVTdGFja1RyYWNlKHguc3RhY2spKSA6ICcnfWA7XHJcbn1cclxuXHJcbmFzeW5jIGZ1bmN0aW9uIGV4ZWNUZXN0KHQ6IFRlc3QsIHByZWRpY2F0ZTogc3RyaW5nIHwgdW5kZWZpbmVkLCBsb2dzOiBhbnlbXSxcclxuICB0ZXN0VGltZW91dD86IG51bWJlciwgcGFja2FnZU5hbWU/OiBzdHJpbmcsIHZlcmJvc2U/OiBib29sZWFuXHJcbik6IFByb21pc2U8VGVzdFJlc3VsdCB8IHVuZGVmaW5lZD4ge1xyXG4gIGxvZ3MubGVuZ3RoID0gMDtcclxuICBsZXQgcjogVGVzdFJlc3VsdDtcclxuICBsZXQgdHlwZTogc3RyaW5nID0gJ3BhY2thZ2UnO1xyXG4gIGNvbnN0IGZpbHRlciA9IHByZWRpY2F0ZSAhPSB1bmRlZmluZWQgJiYgKHQubmFtZS50b0xvd2VyQ2FzZSgpICE9PSBwcmVkaWNhdGUudG9Mb3dlckNhc2UoKSk7XHJcbiAgbGV0IHNraXAgPSB0Lm9wdGlvbnM/LnNraXBSZWFzb24gfHwgZmlsdGVyO1xyXG4gIGxldCBza2lwUmVhc29uID0gZmlsdGVyID8gJ3NraXBwZWQnIDogdC5vcHRpb25zPy5za2lwUmVhc29uO1xyXG5cclxuICBpZiAoREcuVGVzdC5pc0luQmVuY2htYXJrICYmICF0Lm9wdGlvbnM/LmJlbmNobWFyaykge1xyXG4gICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFNraXBwZWQge3ske3QuY2F0ZWdvcnl9fX0ge3ske3QubmFtZX19fSBkb2VzbnQgYXZhaWxhYmxlIGluIGJlbmNobWFyayBtb2RlYCk7XHJcbiAgICByZXR1cm4gdW5kZWZpbmVkO1xyXG4gIH1cclxuXHJcbiAgaWYgKCFza2lwKVxyXG4gICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFN0YXJ0ZWQge3ske3QuY2F0ZWdvcnl9fX0ge3ske3QubmFtZX19fWApO1xyXG4gIGNvbnN0IHN0YXJ0ID0gRGF0ZS5ub3coKTtcclxuICBjb25zdCBzdGFydERhdGUgPSBuZXcgRGF0ZShzdGFydCkudG9JU09TdHJpbmcoKTtcclxuICB0cnkge1xyXG4gICAgaWYgKHNraXApXHJcbiAgICAgIHIgPSB7IG5hbWU6IHQubmFtZSwgb3duZXI6dC5vcHRpb25zPy5vd25lciA/PyAnJywgY2F0ZWdvcnk6IHQuY2F0ZWdvcnksIGxvZ3M6ICcnLCBkYXRlOiBzdGFydERhdGUsIHN1Y2Nlc3M6IHRydWUsIHJlc3VsdDogc2tpcFJlYXNvbiEsIG1zOiAwLCBza2lwcGVkOiB0cnVlLCBwYWNrYWdlOiBwYWNrYWdlTmFtZSA/PyAnJywgZmxha2luZzogREcuVGVzdC5pc1JlcHJvZHVjaW5nfTtcclxuICAgIGVsc2Uge1xyXG4gICAgICBsZXQgdGltZW91dF8gPSB0ZXN0VGltZW91dCA/PyBTVEFOREFSVF9USU1FT1VUO1xyXG5cclxuICAgICAgaWYgKERHLlRlc3QuaXNQcm9maWxpbmcpXHJcbiAgICAgICAgY29uc29sZS5wcm9maWxlKGAke3QuY2F0ZWdvcnl9OiAke3QubmFtZX1gKTtcclxuXHJcbiAgICAgIHIgPSB7IG5hbWU6IHQubmFtZSwgb3duZXI6dC5vcHRpb25zPy5vd25lciA/PyAnJywgY2F0ZWdvcnk6IHQuY2F0ZWdvcnksIGxvZ3M6ICcnLCBkYXRlOiBzdGFydERhdGUsIHN1Y2Nlc3M6IHRydWUsIHJlc3VsdDogKGF3YWl0IHRpbWVvdXQodC50ZXN0LCB0aW1lb3V0XykpLnRvU3RyaW5nKCkgPz8gJ09LJywgbXM6IDAsIHNraXBwZWQ6IGZhbHNlICwgcGFja2FnZTogcGFja2FnZU5hbWUgPz8gJycsIGZsYWtpbmc6IERHLlRlc3QuaXNSZXByb2R1Y2luZ307XHJcblxyXG4gICAgICBpZiAoREcuVGVzdC5pc1Byb2ZpbGluZykge1xyXG4gICAgICAgIGNvbnNvbGUucHJvZmlsZUVuZChgJHt0LmNhdGVnb3J5fTogJHt0Lm5hbWV9YCk7XHJcbiAgICAgICAgZ3Jvay5zaGVsbC5pbmZvKGBQcm9maWxpbmcgb2YgJHt0LmNhdGVnb3J5fTogJHt0Lm5hbWV9IGZpbmlzaGVkIFxcbiBQbGVhc2UgZW5zdXJlIHRoYXQgeW91IGhhdmUgb3BlbmVkIERldlRvb2xzIChGMTIpIC8gUGVyZm9ybWFuY2UgcGFuZWwgYmVmb3JlIHRlc3Qgc3RhcnRzLmApO1xyXG4gICAgICB9XHJcbiAgICB9XHJcbiAgfSBjYXRjaCAoeDogYW55KSB7XHJcbiAgICBzdGRFcnJvcih4KTtcclxuICAgIHIgPSB7IG5hbWU6IHQubmFtZSwgb3duZXI6dC5vcHRpb25zPy5vd25lciA/PyAnJywgY2F0ZWdvcnk6IHQuY2F0ZWdvcnksIGxvZ3M6ICcnLCBkYXRlOiBzdGFydERhdGUsIHN1Y2Nlc3M6IGZhbHNlLCByZXN1bHQ6IGF3YWl0IGdldFJlc3VsdCh4KSwgbXM6IDAsIHNraXBwZWQ6IGZhbHNlLCBwYWNrYWdlOiBwYWNrYWdlTmFtZSA/PyAnJywgZmxha2luZzogZmFsc2V9O1xyXG4gIH1cclxuICBpZiAodC5vcHRpb25zPy5pc0FnZ3JlZ2F0ZWQgJiYgci5yZXN1bHQuY29uc3RydWN0b3IgPT09IERHLkRhdGFGcmFtZSkge1xyXG4gICAgY29uc3QgY29sID0gci5yZXN1bHQuY29sKCdzdWNjZXNzJyk7XHJcbiAgICBpZiAoY29sKVxyXG4gICAgICByLnN1Y2Nlc3MgPSBjb2wuc3RhdHMuc3VtID09PSBjb2wubGVuZ3RoO1xyXG4gICAgaWYgKCF2ZXJib3NlKSB7XHJcbiAgICAgIGNvbnN0IGRmID0gci5yZXN1bHQ7XHJcbiAgICAgIGRmLmNvbHVtbnMucmVtb3ZlKCdzdGFjaycpO1xyXG4gICAgICBkZi5yb3dzLnJlbW92ZVdoZXJlKChyKSA9PiByLmdldCgnc3VjY2VzcycpKTtcclxuICAgICAgci5yZXN1bHQgPSBkZjtcclxuICAgIH1cclxuICAgIHIucmVzdWx0ID0gci5yZXN1bHQudG9Dc3YoKTtcclxuICB9XHJcbiAgci5sb2dzID0gbG9ncy5qb2luKCdcXG4nKTtcclxuICByLm1zID0gRGF0ZS5ub3coKSAtIHN0YXJ0O1xyXG4gIGlmICghc2tpcClcclxuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBGaW5pc2hlZCB7eyR7dC5jYXRlZ29yeX19fSB7eyR7dC5uYW1lfX19IHdpdGgge3ske3Iuc3VjY2VzcyA/ICdzdWNjZXNzJyA6ICdlcnJvcid9fX0gZm9yICR7ci5tc30gbXNgKTtcclxuICBpZiAoIXIuc3VjY2Vzcykge1xyXG4gICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogUmVzdWx0IGZvciB7eyR7dC5jYXRlZ29yeX19fSB7eyR7dC5uYW1lfX19OiAke3IucmVzdWx0fWApO1xyXG4gIH1cclxuICByLmNhdGVnb3J5ID0gdC5jYXRlZ29yeTtcclxuICByLm5hbWUgPSB0Lm5hbWU7XHJcbiAgci5vd25lciA9IHQub3B0aW9ucz8ub3duZXIgPz8gJyc7XHJcbiAgaWYgKCFmaWx0ZXIpIHtcclxuICAgIGxldCBwYXJhbXMgPSB7XHJcbiAgICAgICdzdWNjZXNzJzogci5zdWNjZXNzLCAncmVzdWx0Jzogci5yZXN1bHQsICdtcyc6IHIubXMsICdkYXRlJzogci5kYXRlLFxyXG4gICAgICAnc2tpcHBlZCc6IHIuc2tpcHBlZCwgJ2NhdGVnb3J5JzogdC5jYXRlZ29yeSwgJ25hbWUnOiB0Lm5hbWUsICdsb2dzJzogci5sb2dzLCAnb3duZXInOiByLm93bmVyLFxyXG4gICAgICAnZmxha2luZyc6IERHLlRlc3QuaXNSZXByb2R1Y2luZyAmJiByLnN1Y2Nlc3MsXHJcbiAgICAgICdwYWNrYWdlJzogci5wYWNrYWdlXHJcbiAgICB9O1xyXG4gICAgaWYgKHIucmVzdWx0LmNvbnN0cnVjdG9yID09IE9iamVjdCkge1xyXG4gICAgICBjb25zdCByZXMgPSBPYmplY3Qua2V5cyhyLnJlc3VsdCkucmVkdWNlKChhY2MsIGspID0+ICh7IC4uLmFjYywgWydyZXN1bHQuJyArIGtdOiByLnJlc3VsdFtrXSB9KSwge30pO1xyXG4gICAgICBwYXJhbXMgPSB7IC4uLnBhcmFtcywgLi4ucmVzIH07XHJcbiAgICB9XHJcblxyXG4gICAgaWYgKHBhcmFtcy5yZXN1bHQgaW5zdGFuY2VvZiBERy5EYXRhRnJhbWUpXHJcbiAgICAgIHBhcmFtcy5yZXN1bHQgPSBKU09OLnN0cmluZ2lmeShwYXJhbXMucmVzdWx0Py50b0pzb24oKSkgfHwgJyc7XHJcbiAgICBhd2FpdCBncm9rLnNoZWxsLnJlcG9ydFRlc3QodHlwZSwgcGFyYW1zKTtcclxuICB9XHJcbiAgcmV0dXJuIHI7XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiBzaHVmZmxlKGFycmF5OiBhbnlbXSk6IGFueVtdIHtcclxuICBjb25zdCBuZXdBcnIgPSBhcnJheS5zbGljZSgpO1xyXG4gIG5ld0Fyci5zb3J0KCgpID0+IE1hdGgucmFuZG9tKCkgLSAwLjUpO1xyXG4gIHJldHVybiBuZXdBcnI7XHJcbn1cclxuXHJcbi8qIFdhaXRzIFttc10gbWlsbGlzZWNvbmRzICovXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBkZWxheShtczogbnVtYmVyKSB7XHJcbiAgYXdhaXQgbmV3IFByb21pc2UoKHIpID0+IHNldFRpbWVvdXQociwgbXMpKTtcclxufVxyXG5cclxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGF3YWl0Q2hlY2soY2hlY2tIYW5kbGVyOiAoKSA9PiBib29sZWFuLFxyXG4gIGVycm9yOiBzdHJpbmcgPSAnVGltZW91dCBleGNlZWRlZCcsIHdhaXQ6IG51bWJlciA9IDUwMCwgaW50ZXJ2YWw6IG51bWJlciA9IDUwKTogUHJvbWlzZTxhbnk+IHtcclxuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUsIHJlamVjdCkgPT4ge1xyXG4gICAgc2V0VGltZW91dCgoKSA9PiB7XHJcbiAgICAgIGNsZWFySW50ZXJ2YWwoaW50ZXJ2YWxJZCk7XHJcbiAgICAgIHJlamVjdChuZXcgRXJyb3IoZXJyb3IpKTtcclxuICAgIH0sIHdhaXQpO1xyXG4gICAgLy8gQHRzLWlnbm9yZVxyXG4gICAgY29uc3QgaW50ZXJ2YWxJZDogVGltZW91dCA9IHNldEludGVydmFsKCgpID0+IHtcclxuICAgICAgaWYgKGNoZWNrSGFuZGxlcigpKSB7XHJcbiAgICAgICAgY2xlYXJJbnRlcnZhbChpbnRlcnZhbElkKTtcclxuICAgICAgICByZXNvbHZlKG51bGwpO1xyXG4gICAgICB9XHJcbiAgICB9LCBpbnRlcnZhbCk7XHJcbiAgfSk7XHJcbn1cclxuXHJcbi8vIFJldHVybnMgdGVzdCBleGVjdXRpb24gcmVzdWx0IG9yIGFuIGVycm9yIGluIGNhc2Ugb2YgdGltZW91dFxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdGltZW91dChmdW5jOiAoKSA9PiBQcm9taXNlPGFueT4sIHRlc3RUaW1lb3V0OiBudW1iZXIsIHRpbWVvdXRSZWFzb246IHN0cmluZyA9ICdFWEVDVVRJT04gVElNRU9VVCcpOiBQcm9taXNlPGFueT4ge1xyXG4gIGxldCB0aW1lb3V0OiBhbnkgPSBudWxsO1xyXG4gIGNvbnN0IHRpbWVvdXRQcm9taXNlID0gbmV3IFByb21pc2U8YW55PigoXywgcmVqZWN0KSA9PiB7XHJcbiAgICB0aW1lb3V0ID0gc2V0VGltZW91dCgoKSA9PiB7XHJcbiAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBwcmVmZXItcHJvbWlzZS1yZWplY3QtZXJyb3JzXHJcbiAgICAgIHJlamVjdCh0aW1lb3V0UmVhc29uKTtcclxuICAgIH0sIHRlc3RUaW1lb3V0KTtcclxuICB9KTtcclxuICB0cnkge1xyXG4gICAgcmV0dXJuIGF3YWl0IFByb21pc2UucmFjZShbZnVuYygpLCB0aW1lb3V0UHJvbWlzZV0pO1xyXG4gIH0gZmluYWxseSB7XHJcbiAgICBpZiAodGltZW91dClcclxuICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIGlzRGlhbG9nUHJlc2VudChkaWFsb2dUaXRsZTogc3RyaW5nKTogYm9vbGVhbiB7XHJcbiAgY29uc3QgZGlhbG9ncyA9IERHLkRpYWxvZy5nZXRPcGVuRGlhbG9ncygpO1xyXG4gIGZvciAobGV0IGkgPSAwOyBpIDwgZGlhbG9ncy5sZW5ndGg7IGkrKykge1xyXG4gICAgaWYgKGRpYWxvZ3NbaV0udGl0bGUgPT0gZGlhbG9nVGl0bGUpXHJcbiAgICAgIHJldHVybiB0cnVlO1xyXG4gIH1cclxuICByZXR1cm4gZmFsc2U7XHJcbn1cclxuXHJcbi8qKiBFeHBlY3RzIGFuIGFzeW5jaHJvbm91cyB7QGxpbmsgYWN0aW9ufSB0byB0aHJvdyBhbiBleGNlcHRpb24uIFVzZSB7QGxpbmsgY2hlY2t9IHRvIHBlcmZvcm1cclxuICogZGVlcGVyIGluc3BlY3Rpb24gb2YgdGhlIGV4Y2VwdGlvbiBpZiBuZWNlc3NhcnkuXHJcbiAqIEBwYXJhbSAge2Z1bmN0aW9uKCk6IFByb21pc2U8dm9pZD59IGFjdGlvblxyXG4gKiBAcGFyYW0gIHtmdW5jdGlvbihhbnkpOiBib29sZWFufSBjaGVja1xyXG4gKiBAcmV0dXJuIHtQcm9taXNlPHZvaWQ+fVxyXG4gKi9cclxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGV4cGVjdEV4Y2VwdGlvbkFzeW5jKGFjdGlvbjogKCkgPT4gUHJvbWlzZTx2b2lkPixcclxuICBjaGVjaz86IChleGNlcHRpb246IGFueSkgPT4gYm9vbGVhbik6IFByb21pc2U8dm9pZD4ge1xyXG4gIGxldCBjYXVnaHQ6IGJvb2xlYW4gPSBmYWxzZTtcclxuICBsZXQgY2hlY2tlZDogYm9vbGVhbiA9IGZhbHNlO1xyXG4gIHRyeSB7XHJcbiAgICBhd2FpdCBhY3Rpb24oKTtcclxuICB9IGNhdGNoIChlKSB7XHJcbiAgICBjYXVnaHQgPSB0cnVlO1xyXG4gICAgY2hlY2tlZCA9ICFjaGVjayB8fCBjaGVjayhlKTtcclxuICB9IGZpbmFsbHkge1xyXG4gICAgaWYgKCFjYXVnaHQpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignQW4gZXhjZXB0aW9uIGlzIGV4cGVjdGVkIGJ1dCBub3QgdGhyb3duJyk7XHJcbiAgICBpZiAoIWNoZWNrZWQpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignQW4gZXhwZWN0ZWQgZXhjZXB0aW9uIGlzIHRocm93biwgYnV0IGl0IGRvZXMgbm90IHNhdGlzZnkgdGhlIGNvbmRpdGlvbicpO1xyXG4gIH1cclxufVxyXG5cclxuY29uc3QgY2F0REYgPSBERy5EYXRhRnJhbWUuZnJvbUNvbHVtbnMoW0RHLkNvbHVtbi5mcm9tU3RyaW5ncygnY29sJywgWyd2YWwxJywgJ3ZhbDInLCAndmFsMyddKV0pO1xyXG5cclxuLyoqXHJcbiAqIFVuaXZlcnNhbCB0ZXN0IGZvciB2aWV3ZXJzLiBJdCBzZWFyY2ggdmlld2VycyBpbiBET00gYnkgdGFnczogY2FudmFzLCBzdmcsIGltZywgaW5wdXQsIGgxLCBhXHJcbiAqIEBwYXJhbSAge3N0cmluZ30gdiBWaWV3ZXIgbmFtZVxyXG4gKiBAcGFyYW0gIHtfREcuRGF0YUZyYW1lfSBkZiBEYXRhZnJhbWUgdG8gdXNlLiBTaG91bGQgaGF2ZSBhdCBsZWFzdCAzIHJvd3NcclxuICogQHBhcmFtICB7Ym9vbGVhbn0gb3B0aW9ucy5kZXRlY3RTZW1hbnRpY1R5cGVzIFNwZWNpZnkgd2hldGhlciB0byBkZXRlY3Qgc2VtYW50aWMgdHlwZXMgb3Igbm90XHJcbiAqIEBwYXJhbSAge2Jvb2xlYW59IG9wdGlvbnMucmVhZE9ubHkgSWYgc2V0IHRvIHRydWUsIHRoZSBkYXRhZnJhbWUgd2lsbCBub3QgYmUgbW9kaWZpZWQgZHVyaW5nIHRoZSB0ZXN0XHJcbiAqIEBwYXJhbSAge2Jvb2xlYW59IG9wdGlvbnMuYXJiaXRyYXJ5RGZUZXN0IElmIHNldCB0byBmYWxzZSwgdGVzdCBvbiBhcmJpdHJhcnkgZGF0YWZyYW1lXHJcbiAqIChvbmUgY2F0ZWdvcmljYWwgY29sdW1uKSB3aWxsIG5vdCBiZSBwZXJmb3JtZWRcclxuICogQHBhcmFtICB7b2JqZWN0fSBvcHRpb25zIExpc3Qgb2Ygb3B0aW9ucyAob3B0aW9uYWwpXHJcbiAqIEByZXR1cm4ge1Byb21pc2U8dm9pZD59IFRoZSB0ZXN0IGlzIGNvbnNpZGVyZWQgc3VjY2Vzc2Z1bCBpZiBpdCBjb21wbGV0ZXMgd2l0aG91dCBlcnJvcnNcclxuICovXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiB0ZXN0Vmlld2VyKHY6IHN0cmluZywgZGY6IF9ERy5EYXRhRnJhbWUsIG9wdGlvbnM/OiB7XHJcbiAgZGV0ZWN0U2VtYW50aWNUeXBlcz86IGJvb2xlYW4sIHJlYWRPbmx5PzogYm9vbGVhbiwgYXJiaXRyYXJ5RGZUZXN0PzogYm9vbGVhbixcclxuICBwYWNrYWdlTmFtZT86IHN0cmluZywgYXdhaXRWaWV3ZXI/OiAodmlld2VyOiBfREcuVmlld2VyKSA9PiBQcm9taXNlPHZvaWQ+XHJcbn0pOiBQcm9taXNlPHZvaWQ+IHtcclxuICBjb25zdCBwYWNrYWdlTmFtZSA9IG9wdGlvbnM/LnBhY2thZ2VOYW1lID8/ICcnO1xyXG4gIGlmIChvcHRpb25zPy5kZXRlY3RTZW1hbnRpY1R5cGVzKVxyXG4gICAgYXdhaXQgZ3Jvay5kYXRhLmRldGVjdFNlbWFudGljVHlwZXMoZGYpO1xyXG4gIGNvbnN0IHR2ID0gZ3Jvay5zaGVsbC5hZGRUYWJsZVZpZXcoZGYpO1xyXG5cclxuICB0cnkge1xyXG4gICAgLy8xLiBPcGVuLCBkbyBub3RoaW5nIGFuZCBjbG9zZVxyXG4gICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCk7XHJcbiAgICAvL2luIGNhc2Ugdmlld2VyIHdpdGggYXN5bmMgcmVuZGVyaW5nIC0gd2FpdCBmb3IgcmVuZGVyIHRvIGNvbXBsZXRlXHJcbiAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXHJcbiAgICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQsIHVuZGVmaW5lZCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpO1xyXG5cclxuICAgIC8vMi4gT3BlbiB2aWV3ZXIsIHJ1biBzZWxlY3Rpb24sIGZpbHRlciwgZXRjLiBhbmQgY2xvc2VcclxuICAgIGlmICghb3B0aW9ucz8ucmVhZE9ubHkpIHtcclxuICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgc2VsZWN0RmlsdGVyQ2hhbmdlQ3VycmVudCk7XHJcbiAgICAgIGlmIChvcHRpb25zPy5hd2FpdFZpZXdlcilcclxuICAgICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBzZWxlY3RGaWx0ZXJDaGFuZ2VDdXJyZW50LCBvcHRpb25zIS5hd2FpdFZpZXdlcik7XHJcbiAgICB9XHJcblxyXG4gICAgLy8yLiBPcGVuIHZpZXdlciwgY2hhbmdlIG9wdGlvbnMsIHNhdmUgbGF5b3V0IGFuZCBjbG9zZVxyXG4gICAgbGV0IHByb3BzQW5kTGF5b3V0OiB7IGxheW91dDogYW55LCBzYXZlZFByb3BzOiBhbnkgfSB8IG51bGwgPSBudWxsO1xyXG4gICAgcHJvcHNBbmRMYXlvdXQgPSBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dCk7XHJcbiAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXHJcbiAgICAgIHByb3BzQW5kTGF5b3V0ID0gYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCxcclxuICAgICAgICBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpXHJcblxyXG4gICAgLy8zLiBMb2FkIGxheW91dFxyXG4gICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3TGF5b3V0QXBwbGllZCwgbG9hZExheW91dCwgdW5kZWZpbmVkLCBwcm9wc0FuZExheW91dD8ubGF5b3V0LFxyXG4gICAgICB7IHNhdmVkUHJvcHM6IHByb3BzQW5kTGF5b3V0Py5zYXZlZFByb3BzIH0pO1xyXG4gICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxyXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdMYXlvdXRBcHBsaWVkLCBsb2FkTGF5b3V0LCBvcHRpb25zIS5hd2FpdFZpZXdlcixcclxuICAgICAgICBwcm9wc0FuZExheW91dD8ubGF5b3V0LCB7IHNhdmVkUHJvcHM6IHByb3BzQW5kTGF5b3V0Py5zYXZlZFByb3BzIH0pO1xyXG5cclxuICAgIC8vNC4gT3BlbiB2aWV3ZXIgb24gYXJiaXRhcnkgZGF0YXNldFxyXG4gICAgaWYgKG9wdGlvbnM/LmFyYml0cmFyeURmVGVzdCAhPT0gZmFsc2UpIHtcclxuICAgICAgdHYuZGF0YUZyYW1lID0gY2F0REY7XHJcbiAgICAgIGF3YWl0IGRlbGF5KDUwKTtcclxuICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCk7XHJcbiAgICAgIGlmIChvcHRpb25zPy5hd2FpdFZpZXdlcilcclxuICAgICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCB1bmRlZmluZWQsIG9wdGlvbnMhLmF3YWl0Vmlld2VyKTtcclxuICAgIH1cclxuXHJcbiAgICAvLzUuIENhbGwgcG9zdHBvbmVkIGZpbHRlcmluZ1xyXG4gICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgZmlsdGVyQXN5bmMpO1xyXG4gICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxyXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBmaWx0ZXJBc3luYywgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpO1xyXG5cclxuICB9IGZpbmFsbHkge1xyXG4gICAgLy8gY2xvc2VBbGwoKSBpcyBoYW5kbGluZyBieSBjb21tb24gdGVzdCB3b3JrZmxvd1xyXG4gICAgLy8gZ3Jvay5zaGVsbC5jbG9zZUFsbCgpO1xyXG4gICAgLy8gREcuQmFsbG9vbi5jbG9zZUFsbCgpO1xyXG4gIH1cclxufVxyXG4iXX0=