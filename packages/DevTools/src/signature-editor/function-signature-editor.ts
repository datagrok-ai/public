import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import CodeMirror from 'codemirror';
import {BehaviorSubject} from 'rxjs';
import {
  COMMON_TAG_NAMES,
  DATA_QUERY_VIEW,
  DEFAULT_CATEGORY,
  DIRECTION,
  FUNC_PARAM_FIELDS,
  FUNC_PROPS_FIELD,
  FUNC_PROPS_FIELDS,
  funcParamTypes,
  functionParamsMapping,
  functionPropsCode,
  functionPropsLabels,
  getChoicesByName,
  headerSign,
  helpUrls,
  highlightModeByLang,
  LANGUAGE,
  obligatoryFuncProps,
  OPTIONAL_TAG_NAME,
  OPTIONAL_TAG_NAMES,
  optionTags,
  tooltipMessage
} from './const';
import 'codemirror/mode/javascript/javascript';
import 'codemirror/mode/python/python';
import 'codemirror/mode/octave/octave';
import 'codemirror/mode/r/r';
import 'codemirror/mode/julia/julia';
import 'codemirror/mode/sql/sql.js';
import '../../css/styles.css';

export function functionSignatureEditor(view: DG.View) {
  view.type === DATA_QUERY_VIEW ? addFseRibbonQuery(view) : addFseRibbonScript(view);
}

const editor = (v: DG.View) => (v.root.querySelector('.CodeMirror') as any).CodeMirror.getDoc().getValue();

function addFseRibbonQuery(v: DG.View) {
  setTimeout(() => {
    const iconFse = ui.iconFA('magic', () => openFse(v, editor(v)), 'Open Signature Editor');
    //@ts-ignore
    v.ribbonMenu.root.previousSibling?.append(ui.div(iconFse,'d4-ribbon-item'));
  }, 500);
}

function addFseRibbonScript(v: DG.View) {
  setTimeout(() => {
    const panels = v.getRibbonPanels();
    const iconFse = ui.iconFA('magic', () => openFse(v, editor(v)), 'Open Signature Editor');
    if (!panels.some((panel) => panel.some((icon) => {
      return (icon.firstChild as HTMLElement).outerHTML === iconFse.outerHTML;
    }))) 
    v.setRibbonPanels([...panels, [iconFse]]);
  }, 500);
}

function getInputBaseArray(props: DG.Property[], param: any): DG.InputBase[] {
  return props.map((prop) => {
    const input = DG.InputBase.forProperty(prop, param);
    input.setTooltip(tooltipMessage[prop.name as OPTIONAL_TAG_NAME]);
    input.root.onclick = (_) => {
      grok.shell.windows.showHelp = true;
      grok.shell.windows.help.showHelp(helpUrls[prop.caption]);
    };
    return input;
  });
}

let conn; 

async function getDataQuery(sql: string): Promise<DG.DataQuery> {
  const regex = /--connection:\s*(\S+)/;
  const match = sql.match(regex);
  const connectionName = match ? match[1]: null;
  const connection = await grok.functions.eval(connectionName);
  conn = connectionName;
  return connection.query('sql', sql);
}

async function getInputCode(v: DG.View, functionCode: string) {
  return v.type === DATA_QUERY_VIEW
      ? await getDataQuery(functionCode)
      : DG.Script.create(functionCode);
}

async function openFse(v: DG.View, functionCode: string) {
  const inputScriptCopy = await getInputCode(v, functionCode);
  const language = (inputScriptCopy as DG.Script).language || LANGUAGE.SQL;

  const editorView = DG.View.create();
  editorView.name = v.name;

  const openScript = () => {
    editorView.close();
    grok.shell.addView(v);
    const editor = (v.root.querySelector('.CodeMirror') as any).CodeMirror;
    const doc = editor.getDoc();
    doc.setValue(myCM.getDoc().getValue());
  };

  const updateFuncPropValue = (propName: string, v: any) => {
    (inputScriptCopy as any)[propName] = v;
    refreshPreview();
  };

  const generateParamLine = (param: DG.Property, direction: string) => {
    const optionTagsPreview = optionTags(param)
      .map((tag) => {
        return { tag: tag, val: param.options[tag] };
      })
      .filter((value) => !!value.val)
      .map(({ tag, val }) => `${tag}: ${val}`)
      .concat(...(param.category && param.category !== DEFAULT_CATEGORY) ? [`category: ${param.category}`] : [])
      .join('; ');
    return `${headerSign[language]}${direction}: ${param.propertyType} ${param.name ? `${param.name} ` : ''}${param.defaultValue ? `= ${param.defaultValue} ` : ''}${!!optionTagsPreview.length ? `{${optionTagsPreview}} ` : ''}${param.description ? `[${param.description}]` : ''}\n`;
  };

  const nameToProp = (name, options?: {[key: string]: any}) => {
    const prop = DG.Property.create(name, DG.TYPE.STRING,
    (x: any) => x.options[name],
    (x: any, v) => updateValue(x, name, v), '');
    if (options) prop.fromOptions(options);
    return prop;
  };

  const functionProps: DG.Property[] = FUNC_PROPS_FIELDS.map((name) => nameToProp(name, getChoicesByName[name]));

  const getNewProps = () => {
    const newProps = [];
    for (const prop of Object.values(functionProps)) {
      if (!prop.get(inputScriptCopy) && !(inputScriptCopy as any)[prop.name])
        newProps.push(prop);
    }
    return newProps;
  };

  let addButton = ui.div();
  let inputs: HTMLElement = ui.div();

  const functionPropsInput = (prop: DG.Property) => {
    switch (prop.name) {
      case FUNC_PROPS_FIELD.LANGUAGE:
        return ui.input.choice(
          functionPropsLabels[prop.name as FUNC_PROPS_FIELD],
          {value: prop.get(inputScriptCopy) || (inputScriptCopy as any)[prop.name],
          items: prop.choices});
      // case FUNC_PROPS_FIELDS.TAGS:
      //   return ui.input.multiChoice(
      //     functionPropsLabels(prop.name as FUNC_PROPS_FIELDS),
      //     {value: prop.get(inputScriptCopy) || (inputScriptCopy as any)[prop.name],
      //     items: prop.choices}
      //   );
      default:
        return ui.input.string(functionPropsLabels[prop.name as FUNC_PROPS_FIELD],
          {value: prop.get(inputScriptCopy) || (inputScriptCopy as any)[prop.name]});
    }
  };

  const renderAddPropButton = () => {
    if (!getNewProps().length) return ui.div();

    const onItemClick = async (item: DG.Property) => {
      inputScriptCopy[item.name] = ' ';
      const newInputs = ui.inputs(  
        Object.values(functionProps)
          .filter((prop) => !!prop.get(inputScriptCopy) || !!(inputScriptCopy as any)[prop.name])
          .map((prop) => addFullWidthInput(
            functionPropsInput(prop),
            prop,
          )),
      );
      inputs.replaceWith(newInputs);
      inputs = newInputs;
      refreshPreview();
      const newAddButton = renderAddPropButton();
      addButton.replaceWith(newAddButton);
      addButton = newAddButton;
    };
    const menu = DG.Menu.popup();
    getNewProps().forEach(
      (prop) => menu.item(functionPropsLabels[prop.name as FUNC_PROPS_FIELD], () => onItemClick(prop)),
    );

    const button = ui.button([ui.icons.add(() => { })], () => {
      menu.show();
    });

    const div = ui.div([button]);
    div.style.cssText += 'margin-left: 152px';
    return div;
  };

  const addFullWidthInput = (input: DG.InputBase, prop: DG.Property) => {
    (input.root.lastChild as HTMLElement).style.cssText += 'width: 400px; max-width: inherit;';
    input.onInput.subscribe(() => {
      inputScriptCopy[prop.name] = input.stringValue;
      refreshPreview();
    });
    (input.root as HTMLInputElement).placeholder = 'Enter your value...';
    input.root.appendChild(ui.button(ui.icons.delete(() => { }), () => {
      inputScriptCopy[prop.name] = prop.defaultValue;
      refreshPreview();
      input.root.remove();
      const newAddButton = renderAddPropButton();
      addButton.replaceWith(newAddButton);
      addButton = newAddButton;
    }));
    (input.root.lastChild as HTMLElement).style.cssText += 'display: inline-flex; justify-content: center; flex-direction: column';
    if (obligatoryFuncProps.includes(prop.name)) {
      input.root.style.cssText += 'padding-right: 35.375px';
      (input.root.lastChild as HTMLElement).style.cssText += 'display: none;';
    }
    return input;
  };

  const functionPropsForm = () => {
    addButton = renderAddPropButton();

    inputs = ui.inputs(
      Object.values(functionProps)
        .filter((prop) => !!prop.get(inputScriptCopy) || !!(inputScriptCopy as any)[prop.name])
        .map((prop) => addFullWidthInput(
          functionPropsInput(prop),
          prop,
        )),
    );

    return ui.panel([
      ui.divV([
        inputs,
        addButton,
      ]),
    ]);
  };

  let functionParamsCopy = [
    ...inputScriptCopy.inputs.map((prop) => {
      prop.options.direction = DIRECTION.INPUT; return prop;
    }),
    ...inputScriptCopy.outputs.map((prop) => {
      prop.options.direction = DIRECTION.OUTPUT; return prop;
    }),
  ];

  const onFunctionParamClick = (paramName: string) => {
    const paramIndex = functionParamsCopy.findIndex((input) => input.name === paramName);
    const param = functionParamsCopy[paramIndex];

    if (!param) return ui.div('');

    const obligatoryTagsInputBase = getInputBaseArray(obligatoryFuncParamsTags, param);

    const optionalTagsInputBase = getInputBaseArray(
      optionalFuncParamsTags.filter((prop) => optionTags(param).includes(prop.name as OPTIONAL_TAG_NAME)),
      param
    );

    const result = ui.form([
      ...obligatoryTagsInputBase,
      ...optionalTagsInputBase,
    ]);

    result.style.display = 'flex';
    result.style.flexDirection = 'column';

    const helpIcon = ui.iconFA('question', () => {
      window.open(
        'https://datagrok.ai/help/datagrok/functions/func-params-annotation', 
        '_blank'
      );
    });
    helpIcon.classList.add('dt-help-icon');

    grok.shell.o = ui.divV([
      ui.divH(
        [
          ui.h1(`Param: ${paramName}`), 
          helpIcon
        ],
      ), ui.block75([result])
    ]);
  };

  const updateValue = (param: DG.Property, propName: string, v: any) => {
    (param as any)[propName] = v;
    const globalParam =
      functionParamsCopy.find((copy) => copy.name === param.name);
    if (globalParam) globalParam.options[propName] = v;
    functionParamsState.next(functionParamsCopy);
    refreshPreview();
  };

  const obligatoryFuncParamsProps: DG.Property[] = [
    (() => {
      const temp = DG.Property.create(FUNC_PARAM_FIELDS.DIRECTION, DG.TYPE.STRING,
        (x: any) => x.options?.[FUNC_PARAM_FIELDS.DIRECTION],
        (x: any, v) => updateFuncPropValue(FUNC_PARAM_FIELDS.DIRECTION, v), '');
      temp.fromOptions({ choices: [DIRECTION.INPUT, DIRECTION.OUTPUT] });
      return temp;
    })(),
    DG.Property.create(FUNC_PARAM_FIELDS.NAME, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.NAME],
      (x: any, v) => updateFuncPropValue(FUNC_PARAM_FIELDS.NAME, v), ''),
    (() => {
      const temp = DG.Property.create(FUNC_PARAM_FIELDS.TYPE, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.TYPE],
        (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.TYPE, v), '');
      temp.fromOptions({ choices: funcParamTypes });
      return temp;
    })(),
    DG.Property.create(FUNC_PARAM_FIELDS.DEFAULT_VALUE, DG.TYPE.STRING,
      (x: any) => String(DG.toJs(x)[FUNC_PARAM_FIELDS.DEFAULT_VALUE] || ''),
      (x: any, v) => updateFuncPropValue(FUNC_PARAM_FIELDS.DEFAULT_VALUE, v), ''),
    (() => {
      const temp = DG.Property.create(FUNC_PARAM_FIELDS.DESCRIPTION, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.DESCRIPTION], (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.DESCRIPTION, v), '');
      temp.fromOptions({ editor: 'textarea' });
      return temp;
    })(),
    DG.Property.create(FUNC_PARAM_FIELDS.CATEGORY, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.CATEGORY],
      (x: any, v) => updateFuncPropValue(FUNC_PARAM_FIELDS.CATEGORY, v), ''),
  ];

  const obligatoryFuncParamsTags: DG.Property[] = COMMON_TAG_NAMES.map((name) => nameToProp(name));

  const optionalFuncParamsTags: DG.Property[] = OPTIONAL_TAG_NAMES.map((name) => nameToProp(name, getChoicesByName[name]));
  
  const paramsDF = DG.DataFrame.create(functionParamsCopy.length);
  for (const p of obligatoryFuncParamsProps) {
    (paramsDF.columns as DG.ColumnList)
      .addNew(functionParamsMapping[p.name as keyof typeof functionParamsMapping], p.propertyType as DG.ColumnType)
      .init((i: number) => {
        return p.get(functionParamsCopy[i]) ||
          (functionParamsCopy[i] as any)[p.name] ||
          functionParamsCopy[i].options[p.name];
      });
  }
  (paramsDF.columns as DG.ColumnList).addNew('+', DG.TYPE.STRING);

  paramsDF.onCurrentRowChanged.subscribe(() => {
    if (paramsDF?.currentRow && paramsDF.currentRow.idx !== -1)
      onFunctionParamClick(paramsDF.currentRow.get('Name'));
  });

  const paramsGrid = DG.Grid.create(paramsDF);
  paramsGrid.root.style.width = '100%';
  paramsGrid.dataFrame.getCol(functionParamsMapping[FUNC_PARAM_FIELDS.TYPE as keyof typeof functionParamsMapping])
    .meta.choices = funcParamTypes;
  paramsGrid.dataFrame.getCol(functionParamsMapping[FUNC_PARAM_FIELDS.DIRECTION as keyof typeof functionParamsMapping])
    .meta.choices = [DIRECTION.INPUT, DIRECTION.OUTPUT];
  paramsGrid.setOptions({'showColumnGridlines': false});

  const col = paramsGrid.columns.byName('+');
  col.cellType = 'html';
  col.width = 68;

  paramsGrid.onCellPrepare((gc) => {
    if (gc.gridColumn.name !== '+')
      return;


    const deleteBtn = (rowIdx: number) => ui.div(
      ui.icons.delete(() => {
        functionParamsCopy.splice(rowIdx, 1);
        functionParamsState.next(functionParamsCopy);
      }, 'Remove the param'), { style: { 'text-align': 'center', 'margin': '6px' } },
    );


    if (gc.isTableCell) {
      gc.style.element = ui.divH([
        deleteBtn(gc.gridRow),
        addParamBtn(gc.cell.rowIndex)
      ], {style: {'width': '100px'}});
      gc.style.element.style.display = 'flex';
      gc.style.element.style.justifyContent = 'center';
      gc.style.element.style.color = 'var(--blue-1)';
    }

    if (gc.isColHeader)
      gc.customText = '';
  });

  const codeArea = ui.input.textArea('', {value: ''});
  const myCM = CodeMirror.fromTextArea((codeArea.input as HTMLTextAreaElement), { mode: highlightModeByLang[language] });
  const uiArea = await inputScriptCopy.prepare().getEditor();
  const codePanel = ui.panel([codeArea.root]);
  codePanel.style.height = '100%';
  codeArea.root.style.height = '100%';

  const previewTabs = ui.tabControl(
    {
      'CODE': codePanel,
      'UI': ui.panel([uiArea]),
    }).root;

  previewTabs.style.width = '100%';
  previewTabs.style.flexGrow = '3';

  const propsForm = functionPropsForm();

  const defaultParamName = 'newParam';

  const getUniqueParamName = () => {
    const params = functionParamsCopy.filter((p) => p.name.startsWith(defaultParamName));
    let num = 0;
    for (const param of params) {
      const match = param.name.match(/newParam_(\d+)/);
      if (match || param.name === defaultParamName)
        num++;
    }
    return num === 0 ? defaultParamName : `${defaultParamName}_${num}`
  }

  const addParamBtn = (rowIdx: number) => {
    const btn = ui.button(
      [
        ui.div(ui.icons.add(() => {
        }, 'Add the param'), { style: { 'text-align': 'center', 'margin': '6px'} }),
      ],
      () => {
        const newParam = {
          [FUNC_PARAM_FIELDS.DIRECTION]: DIRECTION.INPUT,
          [FUNC_PARAM_FIELDS.NAME]: getUniqueParamName(),
          [FUNC_PARAM_FIELDS.TYPE]: DG.TYPE.BOOL,
          [FUNC_PARAM_FIELDS.DEFAULT_VALUE]: false,
          [FUNC_PARAM_FIELDS.DESCRIPTION]: '',
          [FUNC_PARAM_FIELDS.CATEGORY]: '',
        };

        const t = DG.Property.create(newParam.name, newParam.propertyType, (x: any) => x, (x: any, v) => x = v);
        t.options.direction = newParam.direction;
        t.description = newParam.description;
        functionParamsCopy.splice(rowIdx + 1, 0, t);
        functionParamsState.next(functionParamsCopy);
      },
    );
    btn.style.margin = 'initial';
    return btn;
  };

  const editorTabs = ui.tabControl({
    'PROPERTIES': propsForm,
    'PARAMETERS': ui.divV([
      paramsGrid.root,
    ]),
  });
  editorTabs.root.style.width = '100%';
  editorTabs.root.style.flexGrow = '3';

  editorView.append(
    ui.divV([
      editorTabs,
      previewTabs,
    ]),
  );
  editorView.box = true;
  editorView.setRibbonPanels([
    [
      ui.iconFA('eye', () => {
        previewTabs.hidden ? previewTabs.hidden = false : previewTabs.hidden = true;
      }),
      ui.iconFA('code', () => openScript(), 'Open function editor'),
    ],
  ]);

  const refreshPreview = async () => {
    let result = '';
    if (v.type === DATA_QUERY_VIEW) 
      result += `${headerSign[language]}${functionPropsCode[FUNC_PROPS_FIELD.CONNECTION]}: ${conn}\n`;
    Object.values(functionProps).map((propField) => {
      const propValue = propField.get(inputScriptCopy) || (inputScriptCopy as any)[propField.name];
      if (!!propValue && !!propValue.length) {
        const propName = functionPropsCode[propField.name];
        result += `${headerSign[language]}${functionPropsCode[propName]}: ${propValue}\n`;
      }
    });
    functionParamsCopy.map((param) => {
      result += generateParamLine(param, param.options.direction);
    });
    const regex = new RegExp(`^(${headerSign[language]}.*\n)*`, 'g');
    result += (v.type === DATA_QUERY_VIEW) 
    ? (inputScriptCopy as DG.DataQuery).query.substring((inputScriptCopy as DG.DataQuery).query.match(regex)[0].length + 1) 
    : (inputScriptCopy as DG.Script).script.substring((inputScriptCopy as DG.Script).script.match(regex)[0].length + 1);
    myCM.setOption('mode', highlightModeByLang[language]);
    myCM.setValue(result);
    myCM.setSize('100%', '100%');

    const modifiedScript = await getInputCode(v, result);
    const newUiArea = await modifiedScript.prepare().getEditor();
    uiArea.innerHTML = '';
    uiArea.append(newUiArea);
  };

  v.close();
  grok.shell.addView(editorView);
  refreshPreview();

  const functionParamsState = new BehaviorSubject(functionParamsCopy);
  functionParamsState.subscribe(() => {
    if (functionParamsCopy.length === paramsDF.rowCount) {
      for (const p of obligatoryFuncParamsProps) {
        (paramsDF.columns as DG.ColumnList)
          .byName(functionParamsMapping[p.name as keyof typeof functionParamsMapping])
          .init((i: number) => {
            return p.get(functionParamsCopy[i]) ||
              (functionParamsCopy[i] as any)[p.name] ||
              functionParamsCopy[i].options[p.name];
          });
      }
    }

    if (functionParamsCopy.length > paramsDF.rowCount) {
      const newProps = functionParamsCopy.map((el: DG.Property) => el.caption);
      const oldProps = paramsDF.columns.byName(FUNC_PARAM_FIELDS.NAME).toList();
      let rowIdx = newProps.map((element, index) => [element, index])
      .filter(([element, index]) => element !== oldProps[index])
      .map(([_, index]) => index);
      const newParam = functionParamsCopy[rowIdx[0]];
      paramsDF.rows.insertAt(rowIdx[0] as number, 1);
      paramsDF.rows.setValues(rowIdx[0] as number, [
        newParam.options.direction, newParam.name,
        newParam.propertyType, newParam.defaultValue,
        newParam.description, newParam.category,
      ],)
    }

    if (functionParamsCopy.length < paramsDF.rowCount) {
      paramsDF.rows.removeAt(0, paramsDF.rowCount);
      functionParamsCopy.forEach((newParam) => {
        paramsDF.rows.addNew(
          [
            newParam.options.direction, newParam.name,
            newParam.propertyType, String(newParam.defaultValue || ''),
            newParam.description, newParam.category,
            '',
          ],
        );
      });
    }
    refreshPreview();
  });

  paramsGrid.onCellValueEdited.subscribe((editedCell) => {
    const rowIndex = functionParamsCopy.findIndex((param) => param.name === (editedCell.tableRow as any)['Name']);
    if (rowIndex) {
      (functionParamsCopy[rowIndex] as any)
      [functionParamsMapping[editedCell.cell.column.name as keyof typeof functionParamsMapping]] = editedCell.cell.value || undefined;
      functionParamsCopy[rowIndex].options
      [functionParamsMapping[editedCell.cell.column.name as keyof typeof functionParamsMapping]] = editedCell.cell.value || undefined;
      functionParamsState.next(functionParamsCopy);
    }
    onFunctionParamClick((editedCell.tableRow as any)['Name']);
  });
}
