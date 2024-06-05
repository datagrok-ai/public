import React, { useEffect } from 'react';

const ToolMarkup = () => {

    const toolMarkup = (<div class="MarkUpTool">
        <label for="compare">Compare</label>
        <select id="selector1" name="compare" style={{margin: '0 5px'}}>
        </select>
        <label for="with"> with </label>
        <select id="selector2" name="with" style={{margin: '0 5px'}}>
        </select>
        <br />
        <input type="checkbox" id="breakingChangesCheckbox" style={{margin: '0 2px'}}/>
        <label for="horns">Breaking changes</label>
        <br />
        <input type="checkbox" id="anyCausesBreakingChangesCheckbox" disabled  style={{margin: '0 2px'}}/>
        <label for="horns">Show type 'any' as breaking changes</label>
        <br />
        <input type="checkbox" id="newFunctionalityCheckbox"  style={{margin: '0 2px'}}/>
        <label for="horns">New functionality</label>
        <br />
        <br />

        <div id="compareResultsBlock">
        </div>
        <script src="./output.js" type="text/javascript"></script> 
    </div>);
    
    useEffect(() => {
        var breakingChangesCheckbox = document.getElementById("breakingChangesCheckbox");


        var selectorToCompare1 = document.getElementById("selector1");
        var selectorToCompare2 = document.getElementById("selector2");
        var compareResultsBlock = document.getElementById("compareResultsBlock");

        var filesUrl = "/versions/";
        var versionsFile = "versions.json";

        var showBreakingChanges = false;
        var showNewFunctionality = false;
        var isAnyTypeCausesBreakingChanges = false;

        var breakingChangesCheckbox = document.getElementById("breakingChangesCheckbox");
        var newFunctionalityCheckbox = document.getElementById("newFunctionalityCheckbox");
        var anyCausesBreakingChangesCheckbox = document.getElementById("anyCausesBreakingChangesCheckbox");

        var promiseRegex = new RegExp("Promise\<([^\'>]\>)");
        var anyRegex = new RegExp("any");

        breakingChangesCheckbox.addEventListener('change', function () {
            if (this.checked) {
                showBreakingChanges = true;
                anyCausesBreakingChangesCheckbox.disabled = false;
            } else {
                showBreakingChanges = false;
                anyCausesBreakingChangesCheckbox.disabled = true;
                anyCausesBreakingChangesCheckbox.checked = false;
                isAnyTypeCausesBreakingChanges = false;
            }
            RenderData();
        });

        anyCausesBreakingChangesCheckbox.addEventListener('change', function () {
            if (this.checked) {
                isAnyTypeCausesBreakingChanges = true;
            } else {
                isAnyTypeCausesBreakingChanges = false;
            }
            RenderData();
        });


        newFunctionalityCheckbox.addEventListener('change', function () {
            if (this.checked) {
                showNewFunctionality = true;
            } else {
                showNewFunctionality = false;
            }
            RenderData();
        });



        var versions = [];
        var loadedVersions = [];
        var functionsData = {};

        async function GetVersion() {
            var jsonResult = await (await fetch(filesUrl + versionsFile)).json();
            (jsonResult.versions).forEach(element => {
                versions[versions.length] = element;
            });
            versions = versions.map(a => a.replace(/\d+/g, n => +n + 100000)).sort()
                .map(a => a.replace(/\d+/g, n => +n - 100000));

        }

        function UpdateSelectors() {
            versions.forEach((element) => {
                var option = document.createElement("option");
                option.value = element;
                option.textContent = element;
                selectorToCompare1.appendChild(option)
                option = document.createElement("option");
                option.value = element;
                option.textContent = element;
                selectorToCompare2.appendChild(option)

                selectorToCompare1.selectedIndex = selectorToCompare1.childElementCount - 1;
                selectorToCompare2.selectedIndex = selectorToCompare2.childElementCount - 2;
            });
        }


        async function LoadVersionsData(versionName) {
            if (!loadedVersions.includes(versionName)) {
                var jsonResult = await (await fetch(filesUrl + versionName + ".json")).json();

                Object.keys(jsonResult).forEach(className => {
                    if (!functionsData.hasOwnProperty(className)) {
                        functionsData[className] = {};
                    }

                    Object.keys(jsonResult[className]).forEach((functionName) => {
                        if (!functionsData[className].hasOwnProperty(functionName)) {
                            functionsData[className][functionName] = {};
                        }
                        functionsData[className][functionName][versionName] = jsonResult[className][functionName];
                    });
                });
                loadedVersions[loadedVersions.length] = versionName;
            }
        }


        selectorToCompare1.addEventListener("change", async function (event) {
            await LoadVersionsData(event.target.value);
            RenderData();
        });

        selectorToCompare2.addEventListener("change", async function (event) {
            await LoadVersionsData(event.target.value);
            RenderData();
        });


        GetVersion().then(async () => {
            UpdateSelectors();
            await LoadVersionsData(selectorToCompare1.value);
            await LoadVersionsData(selectorToCompare2.value);
            RenderData();
        });

        function RenderData() {
            compareResultsBlock.innerHTML = '';

            if (selectorToCompare1.value == selectorToCompare2.value) {
                RenederOneVersion();
            }
            else {
                RenederTwoVersion();
            }
        }


        function RenederOneVersion() {
            Object.keys(functionsData).sort().forEach(className => {
                var detailsBlock = document.createElement("details");
                var summaryBlock = document.createElement("summary");
                summaryBlock.textContent = className;
                detailsBlock.appendChild(summaryBlock)
                Object.keys(functionsData[className]).sort().forEach((functionName) => {
                    if (functionsData[className][functionName].hasOwnProperty(selectorToCompare1.value)) {
                        var preElement1 = document.createElement("pre");
                        preElement1.textContent = functionSignatureToString(functionName, functionsData[className][functionName][selectorToCompare1.value]);

                        detailsBlock.appendChild(preElement1)
                        detailsBlock.appendChild(document.createElement("br"))
                    }
                });
                compareResultsBlock.appendChild(detailsBlock)
            });
        }

        function RenederTwoVersion() {

            Object.keys(functionsData).sort().forEach(className => {
                var detailsBlock = document.createElement("details");
                var summaryBlock = document.createElement("summary");
                var hasCompares = false;
                summaryBlock.textContent = className;
                detailsBlock.appendChild(summaryBlock)

                Object.keys(functionsData[className]).sort().forEach((functionName) => {
                    if (functionsData[className][functionName].hasOwnProperty(selectorToCompare1.value)) {
                        if (functionsData[className][functionName].hasOwnProperty(selectorToCompare2.value)) {
                            var preElement1 = document.createElement("pre");
                            var preElement2 = document.createElement("pre");
                            var textElement1 = functionSignatureToString(functionName, functionsData[className][functionName][selectorToCompare1.value]);
                            var textElement2 = functionSignatureToString(functionName, functionsData[className][functionName][selectorToCompare2.value]);

                            if (textElement1 !== textElement2) {
                                if (showBreakingChanges) {

                                    if (functionName.trim() == "forEntity")
                                        debugger
                                    if (isBreakingChanges(functionsData[className][functionName][selectorToCompare1.value],
                                        functionsData[className][functionName][selectorToCompare2.value])) {
                                        detailsBlock.appendChild(preElement1)
                                        detailsBlock.appendChild(preElement2)
                                        detailsBlock.appendChild(document.createElement("br"))
                                        hasCompares = true;
                                    }
                                }
                                else if (showNewFunctionality == false) {
                                    detailsBlock.appendChild(preElement1)
                                    detailsBlock.appendChild(preElement2)
                                    detailsBlock.appendChild(document.createElement("br"))
                                    hasCompares = true;
                                }
                                else {

                                }
                            }
                            preElement1.textContent = selectorToCompare1.value + " : " + textElement1;
                            preElement2.textContent = selectorToCompare2.value + " : " + textElement2;
                        }
                        else {
                            var label = document.createElement("label");
                            label.style = "color: green";
                            label.textContent = "Added";
                            var preElement1 = document.createElement("pre");
                            preElement1.textContent = functionSignatureToString(functionName, functionsData[className][functionName][selectorToCompare1.value]);

                            if (showNewFunctionality) {
                                detailsBlock.appendChild(label)
                                detailsBlock.appendChild(preElement1)
                                detailsBlock.appendChild(document.createElement("br"))
                                hasCompares = true;
                            }
                            else if (showBreakingChanges == false) {
                                detailsBlock.appendChild(label)
                                detailsBlock.appendChild(preElement1)
                                detailsBlock.appendChild(document.createElement("br"))
                                hasCompares = true;
                            }
                        }
                    }
                    else if (functionsData[className][functionName].hasOwnProperty(selectorToCompare2.value)) {
                        var label = document.createElement("label");
                        label.style = "color: red";
                        label.textContent = "Removed";
                        var preElement1 = document.createElement("pre");
                        preElement1.textContent = functionSignatureToString(functionName, functionsData[className][functionName][selectorToCompare2.value]);

                        if (showBreakingChanges) {
                            detailsBlock.appendChild(label)
                            detailsBlock.appendChild(preElement1)
                            detailsBlock.appendChild(document.createElement("br"))
                            hasCompares = true;
                        }
                        else if (showNewFunctionality == false) {
                            detailsBlock.appendChild(label)
                            detailsBlock.appendChild(preElement1)
                            detailsBlock.appendChild(document.createElement("br"))
                            hasCompares = true;
                        }
                    }
                });

                if (hasCompares) {
                    compareResultsBlock.appendChild(detailsBlock)
                }
            });
        }

        function functionSignatureToString(name, data) {
            var result = name + "(";

            if ((data.hasOwnProperty("type")) && data["type"].trim() !== 'get' && data["type"].trim() !== 'set') {
                result = data["type"] + " " + result
            }

            if ((data.hasOwnProperty("async") && data["async"] === true)) {
                result = "async " + result
            }

            if ((data.hasOwnProperty("static") && data["static"] === true)) {
                result = "static " + result
            }

            if ((data.hasOwnProperty("params"))) {
                result = result + data["params"]
            }
            result = result + ")";
            if ((data.hasOwnProperty("result"))) {
                result = result + ": " + data["result"]
            }

            return result;
        }

        function isBreakingChanges(newVersion, oldVersion) {
            if (newVersion["async"] !== oldVersion["async"]) {
                return true;
            }
            if (newVersion["static"] !== oldVersion["static"]) {
                return true;
            }
            if (newVersion["result"] !== oldVersion["result"] && isTypesHasBreakingChanges(newVersion["result"], oldVersion["result"])) {
                return true;
            }
            var oldParams = {};
            var newParams = {};

            if ((oldVersion.hasOwnProperty("params"))) {
                oldVersion["params"].split(",").forEach((element) => {
                    var paramsSplited = element.trim().split(":");
                    oldParams[paramsSplited[0].trim()] = (paramsSplited[1] || '').trim();
                });
            }
            if ((newVersion.hasOwnProperty("params"))) {
                newVersion["params"].split(",").forEach((element) => {
                    var paramsSplited = element.trim().split(":");
                    newParams[paramsSplited[0].trim()] = (paramsSplited[1] || '').trim();
                });
            }

            if (oldParams.length > newParams.length) {
                return true;
            }

            const oldKeys = Object.keys(oldParams);
            const newKeys = Object.keys(newParams);
            let outOfAllFunctions = false;
            for (let i = 0; i < newKeys.length; i++) {
                if (oldKeys.length <= i)
                    outOfAllFunctions = true;

                if (outOfAllFunctions) {
                    if (!newKeys[i].includes("?"))
                        return true;
                }
                else {
                    let isBreakingChangesSignature = isTypesHasBreakingChanges(newParams[newKeys[i]], oldParams[oldKeys[i]]);
                    if (isBreakingChangesSignature)
                        return true;
                }
            }
            return false;
        }

        function isTypesHasBreakingChanges(newType, oldType) {

            var result = false;
            if (newType === oldType)
                return result;

            if (oldType === 'void')
                return result;
            let newHasPromise = promiseRegex.test(newType);
            let oldHasPromise = promiseRegex.test(oldType);


            if ((newHasPromise && !oldHasPromise) || (!newHasPromise && oldHasPromise))
                return true;

            if (!newHasPromise)
                return isTypesSemanticBreakingChanges(newType, oldType);
            return result;
        }


        function isTypesSemanticBreakingChanges(newType, oldType) {
            var result = false;
            if (newType !== oldType) {
                anyRegex.lastIndex = 0;
                let newHasAny = anyRegex.test(newType);
                anyRegex.lastIndex = 0;
                let oldHasAny = anyRegex.test(oldType)

                {
                    let splittedNew = newType.split('|');
                    let splittedOld = oldType.split('|');
                    if (splittedOld.length > splittedNew.length)
                        return true;
                    for (let i = 0; i < splittedNew.length; i++) {
                        splittedNew[i] = splittedNew[i].trim();
                        if (splittedNew[i] === "any") {
                            return isAnyTypeCausesBreakingChanges;
                        }
                    }
                    for (let i = 0; i < splittedOld.length; i++) {
                        if (splittedOld[i] === "any") {
                            return isAnyTypeCausesBreakingChanges;
                        }
                        if (!splittedNew.includes(splittedOld[i].trim())) {
                            return true;
                        }
                    }
                }
            }
            return result;
        }
    });
 
    return toolMarkup;
}

export default ToolMarkup