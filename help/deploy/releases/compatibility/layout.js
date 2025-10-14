import React, { useEffect } from "react";

const ToolMarkup = () => {
  const toolMarkup = (
    <div className="MarkUpTool">
      <label htmlFor="compare">Compare</label>
      <select
        id="selector1"
        name="compare"
        style={{ margin: "0 5px" }}
      ></select>
      <label htmlFor="with"> with </label>
      <select id="selector2" name="with" style={{ margin: "0 5px" }}></select>
      <div
        id="copyButton"
        style={{
          display: "inline-block",
          width: "15px",
          position: "relative",
          marginLeft: "10px",
        }}
      >
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 448 512">
          <path
            fill="#2083d5"
            d="M384 336l-192 0c-8.8 0-16-7.2-16-16l0-256c0-8.8 7.2-16 16-16l140.1 0L400 115.9 400 320c0 8.8-7.2 16-16 16zM192 384l192 0c35.3 0 64-28.7 64-64l0-204.1c0-12.7-5.1-24.9-14.1-33.9L366.1 14.1c-9-9-21.2-14.1-33.9-14.1L192 0c-35.3 0-64 28.7-64 64l0 256c0 35.3 28.7 64 64 64zM64 128c-35.3 0-64 28.7-64 64L0 448c0 35.3 28.7 64 64 64l192 0c35.3 0 64-28.7 64-64l0-32-48 0 0 32c0 8.8-7.2 16-16 16L64 464c-8.8 0-16-7.2-16-16l0-256c0-8.8 7.2-16 16-16l32 0 0-48-32 0z"
          />
        </svg>
      </div>
      <br />
      <input
        type="checkbox"
        id="breakingChangesCheckbox"
        style={{ margin: "0 2px" }}
      />
      <label htmlFor="horns">Breaking changes</label>
      <br />
      <input
        type="checkbox"
        id="anyCausesBreakingChangesCheckbox"
        style={{ margin: "0 2px" }}
      />
      <label htmlFor="horns" style={{}}>
        Type conversions{" "}
      </label>
      <div
        id="tooltip"
        style={{
          display: "inline-block",
          width: "15px",
          position: "relative",
          marginLeft: "10px",
        }}
      >
        <svg id="svgIcon" viewBox="0 0 14 16" style={{ width: "14px" }}>
          <path
            fillRule="evenodd"
            d="M6.3 5.69a.942.942 0 0 1-.28-.7c0-.28.09-.52.28-.7.19-.18.42-.28.7-.28.28 0 .52.09.7.28.18.19.28.42.28.7 0 .28-.09.52-.28.7a1 1 0 0 1-.7.3c-.28 0-.52-.11-.7-.3zM8 7.99c-.02-.25-.11-.48-.31-.69-.2-.19-.42-.3-.69-.31H6c-.27.02-.48.13-.69.31-.2.2-.3.44-.31.69h1v3c.02.27.11.5.31.69.2.2.42.31.69.31h1c.27 0 .48-.11.69-.31.2-.19.3-.42.31-.69H8V7.98v.01zM7 2.3c-3.14 0-5.7 2.54-5.7 5.68 0 3.14 2.56 5.7 5.7 5.7s5.7-2.55 5.7-5.7c0-3.15-2.56-5.69-5.7-5.69v.01zM7 .98c3.86 0 7 3.14 7 7s-3.14 7-7 7-7-3.12-7-7 3.14-7 7-7z"
          ></path>
        </svg>
      </div>
      <br />
      <input
        type="checkbox"
        id="newFunctionalityCheckbox"
        style={{ margin: "0 2px" }}
      />
      <label htmlFor="horns">New functionality</label>
      <br />
      <br />
      <p>New features:</p>
      <ul>
        <li id="classesEncounter">Classes: </li>
        <li id="functionsEncounter">Functions: </li>
      </ul>
      <br />
      <button id="expandAllBtn">Expand All</button>
      <br />

      <div id="compareResultsBlock"></div>
      <script src="./output.js" type="text/javascript"></script>
    </div>
  );

  useEffect(() => {
    if (typeof window !== undefined) {
      const queryParams = new URLSearchParams(window.location.search);

      const from = queryParams.get("from");
      const to = queryParams.get("to");
      const breakingChages = queryParams.get("breakingchanges") ?? "";
      const typeConversions = queryParams.get("typeconversions") ?? "";
      const newFunctionality = queryParams.get("newfunctionality") ?? "";

      var tooltipTrigger = document.getElementById("tooltip");
      var svgIcon = document.getElementById("svgIcon");
      var functionsEncounter = document.getElementById("functionsEncounter");
      var classesEncounter = document.getElementById("classesEncounter");
      var mouseOnSvg = false;
      var expandBtn = document.getElementById("expandAllBtn");
      expandBtn.style = "border-radius: 5px; position: absolute; right: 0;";
      expandBtn.addEventListener("click", () => {
        const details = document.querySelectorAll("details");
        details.forEach((detail) => {
          detail.open = true;
        });
      });

      tooltipTrigger.addEventListener("mouseover", function () {
        var tooltip = document.createElement("span");
        tooltip.textContent =
          "Include cases where parameter type changes either to or from ‘any’";

        tooltip.style.visibility = "visible";
        tooltip.style.width = "520px";
        tooltip.style.backgroundColor = "#555";
        tooltip.style.color = "#fff";
        tooltip.style.textAlign = "center";
        tooltip.style.borderRadius = "6px";
        tooltip.style.padding = "5px 0";
        tooltip.style.position = "absolute";
        tooltip.style.zIndex = "1";
        tooltip.style.bottom = "-10%";
        tooltip.style.left = "70px";
        tooltip.style.marginLeft = "-60px";
        tooltip.style.opacity = "0.9";
        tooltip.style.padding = "5px";
        tooltip.style.transition = "opacity 0.3s";

        tooltipTrigger.appendChild(tooltip);
      });

      var copyButtonTrigger = document.getElementById("copyButton");
      copyButtonTrigger.addEventListener("mouseover", function () {
        var tooltip = document.createElement("span");
        tooltip.textContent = "Copy link";

        tooltip.style.visibility = "visible";
        tooltip.style.width = "120px";
        tooltip.style.backgroundColor = "#555";
        tooltip.style.color = "#fff";
        tooltip.style.textAlign = "center";
        tooltip.style.borderRadius = "6px";
        tooltip.style.padding = "5px 0";
        tooltip.style.position = "absolute";
        tooltip.style.zIndex = "1";
        tooltip.style.bottom = "-10%";
        tooltip.style.left = "20px";
        tooltip.style.marginLeft = "10px";
        tooltip.style.opacity = "0.9";
        tooltip.style.padding = "5px";
        tooltip.style.transition = "opacity 0.3s";

        copyButtonTrigger.appendChild(tooltip);
      });

      tooltipTrigger.addEventListener("mouseout", function () {
        var tooltip = tooltipTrigger.querySelector("span");
        if (tooltip && !mouseOnSvg) {
          tooltip.remove();
        }
      });

      copyButton.addEventListener("mouseout", function () {
        var tooltip = copyButtonTrigger.querySelector("span");
        if (tooltip && !mouseOnSvg) {
          tooltip.remove();
        }
      });

      copyButton.addEventListener("click", function () {
        copyLink();
        copyButton.children[0].children[0].setAttribute("fillRule", "#50A9C5");
        setTimeout(() => {
          copyButton.children[0].children[0].setAttribute(
            "fillRule",
            "#50A9C5"
          );
        }, 500);
      });

      svgIcon.addEventListener("mouseover", function (event) {
        mouseOnSvg = true;
      });

      svgIcon.addEventListener("mouseout", function () {
        mouseOnSvg = false;
      });

      var breakingChangesCheckbox = document.getElementById(
        "breakingChangesCheckbox"
      );

      var selectorToCompare1 = document.getElementById("selector1");
      var selectorToCompare2 = document.getElementById("selector2");
      var compareResultsBlock = document.getElementById("compareResultsBlock");

      var filesUrl = "/versions/";
      var versionsFile = "versions.json";

      var breakingChangesCheckbox = document.getElementById(
        "breakingChangesCheckbox"
      );
      var newFunctionalityCheckbox = document.getElementById(
        "newFunctionalityCheckbox"
      );
      var typeConversionCheckbox = document.getElementById(
        "anyCausesBreakingChangesCheckbox"
      );
      if (typeConversions === "true") typeConversionCheckbox.checked = true;
      if (breakingChages === "true") breakingChangesCheckbox.checked = true;
      if (newFunctionality === "true") newFunctionalityCheckbox.checked = true;

      var promiseRegex = new RegExp("Promise<([^'>]>)");
      var anyRegex = new RegExp("any");

      breakingChangesCheckbox.addEventListener("change", function () {
        RenderData();
      });

      typeConversionCheckbox.addEventListener("change", function () {
        RenderData();
      });

      newFunctionalityCheckbox.addEventListener("change", function () {
        RenderData();
      });

      var versions = [];
      var loadedVersions = [];
      var functionsData = {};
      var addedFunctionsCount = 0;
      var addedClasses = new Set();

      async function GetVersion() {
        var jsonResult = await (await fetch(filesUrl + versionsFile)).json();
        jsonResult.versions.forEach((element) => {
          versions[versions.length] = element;
        });
        versions = versions
          .map((a) => a.replace(/\d+/g, (n) => +n + 100000))
          .sort()
          .map((a) => a.replace(/\d+/g, (n) => +n - 100000));
      }

      function UpdateSelectors() {
        versions.forEach((element) => {
          var option = document.createElement("option");
          option.value = element;
          option.textContent = element;
          selectorToCompare1.appendChild(option);
          option = document.createElement("option");
          option.value = element;
          option.textContent = element;
          selectorToCompare2.appendChild(option);

          selectorToCompare1.selectedIndex =
            selectorToCompare1.childElementCount - 1;
          selectorToCompare2.selectedIndex =
            selectorToCompare2.childElementCount - 2;
        });
      }

      async function LoadVersionsData(versionName) {
        if (!loadedVersions.includes(versionName)) {
          var jsonResult = await (
            await fetch(filesUrl + versionName + ".json")
          ).json();

          Object.keys(jsonResult).forEach((classNameSource) => {
            var className = classNameSource.trim();
            if (!functionsData.hasOwnProperty(className)) {
              functionsData[className] = {};
              functionsData[className].versions = new Set();
            }
            functionsData[className].versions.add(versionName);
            Object.keys(jsonResult[className]).forEach((functionName) => {
              if (!functionsData[className].hasOwnProperty(functionName)) {
                functionsData[className][functionName] = {};
              }
              functionsData[className][functionName][versionName] =
                jsonResult[className][functionName];
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
        if (to && versions.includes(to)) selectorToCompare1.value = to;
        if (from && versions.includes(from)) selectorToCompare2.value = from;
        await LoadVersionsData(selectorToCompare1.value);
        await LoadVersionsData(selectorToCompare2.value);
        RenderData();
      });

      function composeLink() {
        const baseUrl = `${window.location.protocol}//${window.location.host}`;
        const fullPath = `${baseUrl}${window.location.pathname}`;
        let composedLink = `${fullPath}?from=${encodeURIComponent(
          selectorToCompare2.value
        )}&to=${encodeURIComponent(selectorToCompare1.value)}`;

        if (breakingChangesCheckbox.checked)
          composedLink = `${composedLink}&breakingchanges=true`;
        if (typeConversionCheckbox.checked)
          composedLink = `${composedLink}&typeconversions=true`;
        if (newFunctionalityCheckbox.checked)
          composedLink = `${composedLink}&newfunctionality=true`;
        return composedLink;
      }

      function copyLink() {
        var composedLink = composeLink();
        navigator.clipboard
          .writeText(composedLink)
          .then(() => {
            console.log("Composed Link:", composedLink);
          })
          .catch((err) => {
            console.error("Error copying link: ", err);
          });
      }

      function RenderData() {
        compareResultsBlock.innerHTML = "";

        if (selectorToCompare1.value == selectorToCompare2.value) {
          RenederOneVersion();
          classesEncounter.innerText = "Classes: 0";
          functionsEncounter.innerText = "Functions: 0";
        } else {
          RenederTwoVersion();
        }

        history.replaceState(null, "", composeLink());
      }

      function RenederOneVersion() {
        addedClasses = new Set();
        addedFunctionsCount = 0;
        Object.keys(functionsData)
          .sort()
          .forEach((className) => {
            var detailsBlock = document.createElement("details");
            var summaryBlock = document.createElement("summary");
            var data = className.split(" ");
            console.log(className);
            summaryBlock.textContent = data[1];
            detailsBlock.appendChild(summaryBlock);
            Object.keys(functionsData[className])
              .sort()
              .forEach((functionName) => {
                if (
                  functionsData[className][functionName].hasOwnProperty(
                    selectorToCompare1.value
                  ) &&
                  functionsData[className].versions.has(
                    selectorToCompare1.value
                  )
                ) {
                  var preElement1 = document.createElement("pre");
                  preElement1.textContent = functionSignatureToString(
                    functionName,
                    functionsData[className][functionName][
                      selectorToCompare1.value
                    ]
                  );

                  detailsBlock.appendChild(preElement1);
                  detailsBlock.appendChild(document.createElement("br"));
                }
              });
            compareResultsBlock.appendChild(detailsBlock);
          });
      }

      function RenederTwoVersion() {
        addedClasses = new Set();
        addedFunctionsCount = 0;
        Object.keys(functionsData)
          .sort()
          .forEach((className) => {
            var detailsBlock = document.createElement("details");
            var summaryBlock = document.createElement("summary");
            var hasCompares = false;
            var data = className.split(" ");
            summaryBlock.textContent = data[1];
            detailsBlock.appendChild(summaryBlock);
            var addedElements = [];
            var removedElements = [];

            if (
              data[0] === "class" &&
              !functionsData[className].versions.has(
                selectorToCompare2.value
              ) &&
              functionsData[className].versions.has(selectorToCompare1.value)
            ) {
              addedClasses.add(className);
              summaryBlock.style = "color: green";
              console.log(className);
            }
            Object.keys(functionsData[className])
              .sort()
              .forEach((functionName) => {
                if (
                  functionsData[className][functionName].hasOwnProperty(
                    selectorToCompare1.value
                  )
                ) {
                  if (
                    functionsData[className][functionName].hasOwnProperty(
                      selectorToCompare2.value
                    )
                  ) {
                    var preElement1 = document.createElement("pre");
                    var textElement1 = functionSignatureToString(
                      functionName,
                      functionsData[className][functionName][
                        selectorToCompare1.value
                      ]
                    );
                    var textElement2 = functionSignatureToString(
                      functionName,
                      functionsData[className][functionName][
                        selectorToCompare2.value
                      ]
                    );

                    if (textElement1 !== textElement2) {
                      if (
                        breakingChangesCheckbox.checked ||
                        typeConversionCheckbox.checked
                      ) {
                        if (functionName.trim() == "forEntity") debugger;
                        if (
                          isBreakingChanges(
                            functionsData[className][functionName][
                              selectorToCompare1.value
                            ],
                            functionsData[className][functionName][
                              selectorToCompare2.value
                            ]
                          )
                        ) {
                          detailsBlock.appendChild(preElement1);
                          detailsBlock.appendChild(
                            document.createElement("br")
                          );
                          hasCompares = true;
                        }
                      } else if (newFunctionalityCheckbox.checked == false) {
                        detailsBlock.appendChild(preElement1);
                        detailsBlock.appendChild(document.createElement("br"));
                        hasCompares = true;
                      } else {
                      }
                    }

                    var divVersionElement = document.createElement("div");
                    divVersionElement.textContent =
                      selectorToCompare1.value + ": ";
                    divVersionElement.style =
                      "color: green; width: 70px; display: inline-block;";
                    var divSignatureElement = document.createElement("div");
                    divSignatureElement.style = "display: inline-block;";
                    divSignatureElement.textContent = textElement1;

                    preElement1.textContent = "";
                    preElement1.appendChild(divVersionElement);
                    preElement1.appendChild(divSignatureElement);

                    divVersionElement = document.createElement("div");
                    divVersionElement.textContent =
                      selectorToCompare2.value + ": ";
                    divVersionElement.style =
                      "color: red; width: 70px; display: inline-block;";
                    divSignatureElement = document.createElement("div");
                    divSignatureElement.style = "display: inline-block;";
                    divSignatureElement.textContent = textElement2;

                    preElement1.appendChild(document.createElement("br"));
                    preElement1.appendChild(divVersionElement);
                    preElement1.appendChild(divSignatureElement);
                  } else {
                    var preElement1 = document.createElement("pre");
                    preElement1.textContent = functionSignatureToString(
                      functionName,
                      functionsData[className][functionName][
                        selectorToCompare1.value
                      ]
                    );

                    if (newFunctionalityCheckbox.checked) {
                      addedElements[addedElements.length] = preElement1;
                      hasCompares = true;
                    } else if (
                      breakingChangesCheckbox.checked == false &&
                      typeConversionCheckbox.checked == false
                    ) {
                      addedElements[addedElements.length] = preElement1;
                      hasCompares = true;
                    }
                  }
                } else if (
                  functionsData[className][functionName].hasOwnProperty(
                    selectorToCompare2.value
                  )
                ) {
                  var preElement1 = document.createElement("pre");
                  preElement1.textContent = functionSignatureToString(
                    functionName,
                    functionsData[className][functionName][
                      selectorToCompare2.value
                    ]
                  );

                  if (breakingChangesCheckbox.checked) {
                    removedElements[removedElements.length] = preElement1;
                    hasCompares = true;
                  } else if (
                    newFunctionalityCheckbox.checked == false &&
                    typeConversionCheckbox.checked == false
                  ) {
                    removedElements[removedElements.length] = preElement1;
                    hasCompares = true;
                  }
                }
              });

            if (addedElements.length > 0) {
              var table = document.createElement("table");
              var tr = document.createElement("tr");
              var td = document.createElement("td");
              var td2 = document.createElement("td");
              var label = document.createElement("label");

              table.appendChild(tr);
              tr.appendChild(td);
              tr.appendChild(td2);

              table.style = "width: 100%";
              tr.style = "vertical-align: top; border: 0px;";
              td.style = "vertical-align: top; border: 0px;";
              td2.style = "vertical-align: top; border: 0px; width: 100%";

              label.style = "color: green";
              label.textContent = "Added";
              addedFunctionsCount++;
              detailsBlock.appendChild(table);
              td.appendChild(label);
              // detailsBlock.appendChild(document.createElement("br"));

              addedElements.forEach((preElement1) => {
                td2.appendChild(preElement1);
              });
            }

            if (removedElements.length > 0) {
              var label = document.createElement("label");
              label.style = "color: red";
              label.textContent = "Removed";

              detailsBlock.appendChild(label);
              removedElements.forEach((preElement1) => {
                detailsBlock.appendChild(preElement1);
              });
            }

            if (hasCompares) {
              compareResultsBlock.appendChild(detailsBlock);
            }
          });
        classesEncounter.innerText = `Classes: ${addedClasses.size ?? 0}`;
        functionsEncounter.innerText = `Functions: ${addedFunctionsCount}`;
        addedClasses.clear();
      }

      function functionSignatureToString(name, data) {
        var result = name + "(";

        if (
          data.hasOwnProperty("type") &&
          data["type"].trim() !== "get" &&
          data["type"].trim() !== "set"
        ) {
          result = data["type"] + " " + result;
        }

        if (data.hasOwnProperty("async") && data["async"] === true) {
          result = "async " + result;
        }

        if (data.hasOwnProperty("static") && data["static"] === true) {
          result = "static " + result;
        }

        if (data.hasOwnProperty("params")) {
          result = result + data["params"];
        }
        result = result + ")";
        if (data.hasOwnProperty("result")) {
          result = result + ": " + data["result"];
        }

        return result.trim();
      }

      function isBreakingChanges(newVersion, oldVersion) {
        if (newVersion["async"] !== oldVersion["async"]) {
          return breakingChangesCheckbox.checked;
        }
        if (newVersion["static"] !== oldVersion["static"]) {
          return breakingChangesCheckbox.checked;
        }
        if (
          newVersion["result"] !== oldVersion["result"] &&
          isTypesHasBreakingChanges(newVersion["result"], oldVersion["result"])
        ) {
          return breakingChangesCheckbox.checked;
        }
        var oldParams = {};
        var newParams = {};
        if (oldVersion.hasOwnProperty("params")) {
          oldVersion["params"].split(",").forEach((element) => {
            var paramsSplited = element.trim().split(":");
            oldParams[paramsSplited[0].trim()] = (
              paramsSplited[1] || ""
            ).trim();
          });
        }
        if (newVersion.hasOwnProperty("params")) {
          newVersion["params"].split(",").forEach((element) => {
            var paramsSplited = element.trim().split(":");
            newParams[paramsSplited[0].trim()] = (
              paramsSplited[1] || ""
            ).trim();
          });
        }

        if (oldParams.length > newParams.length) {
          return breakingChangesCheckbox.checked;
        }

        const oldKeys = Object.keys(oldParams).filter((str) => str !== "");
        const newKeys = Object.keys(newParams).filter((str) => str !== "");
        let outOfAllFunctions = false;
        for (let i = 0; i < newKeys.length; i++) {
          if (oldKeys.length <= i) outOfAllFunctions = true;

          if (outOfAllFunctions) {
            if (!newKeys[i].includes("?"))
              return breakingChangesCheckbox.checked;
          } else {
            let isBreakingChangesSignature = isTypesHasBreakingChanges(
              newParams[newKeys[i]],
              oldParams[oldKeys[i]]
            );
            if (isBreakingChangesSignature) return true;
          }
        }
        return false;
      }

      function isTypesHasBreakingChanges(newType, oldType) {
        var result = false;
        if (newType === oldType)
          return breakingChangesCheckbox.checked && result;

        if (oldType === "void")
          return breakingChangesCheckbox.checked && result;
        let newHasPromise = promiseRegex.test(newType);
        let oldHasPromise = promiseRegex.test(oldType);

        if (
          (newHasPromise && !oldHasPromise) ||
          (!newHasPromise && oldHasPromise)
        )
          return breakingChangesCheckbox.checked && true;

        if (!newHasPromise)
          return isTypesSemanticBreakingChanges(newType, oldType);
        return breakingChangesCheckbox.checked && result;
      }

      function isTypesSemanticBreakingChanges(newType, oldType) {
        var result = false;
        if (newType !== oldType) {
          anyRegex.lastIndex = 0;
          let newHasAny = anyRegex.test(newType);
          anyRegex.lastIndex = 0;
          let oldHasAny = anyRegex.test(oldType);

          {
            let splittedNew = newType.split("|");
            let splittedOld = oldType.split("|");
            for (let i = 0; i < splittedNew.length; i++) {
              splittedNew[i] = splittedNew[i].trim();
              if (splittedNew[i] === "any") {
                return typeConversionCheckbox.checked;
              }
            }
            for (let i = 0; i < splittedOld.length; i++) {
              if (splittedOld[i] === "any") {
                return typeConversionCheckbox.checked;
              }
              if (!splittedNew.includes(splittedOld[i].trim())) {
                return breakingChangesCheckbox.checked && true;
              }
            }
            if (splittedOld.length > splittedNew.length)
              return breakingChangesCheckbox.checked && true;
          }
        }
        return breakingChangesCheckbox.checked && result;
      }
      RenderData();
    }
  });

  return toolMarkup;
};

export default ToolMarkup;
