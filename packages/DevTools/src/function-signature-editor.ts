import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import CodeMirror from 'codemirror';
import {BehaviorSubject} from 'rxjs';

export async function functionSignatureEditor(view: DG.View) {
    addFseRibbon(view);
}

function addFseRibbon(v: DG.View){
    setTimeout(() => {
      const panels = v.getRibbonPanels();
      // @ts-ignore
      const iconFse = ui.iconFA('magic', () => openFse(v, v.root.lastChild.lastChild.CodeMirror.getDoc().getValue()),'Open Signature Editor');
      if (!panels.some((panel) => panel.some((icon) => {
        return (icon.firstChild as HTMLElement).outerHTML === iconFse.outerHTML
      }))) 
      
      v.setRibbonPanels([...panels, [iconFse]])
    }, 500)
}

const DEFAULT_CATEGORY = 'Misc';

const enum FUNC_PROPS_FIELDS {
  NAME = 'name',
  DESCRIPTION = 'description',
  LANGUAGE = 'language',
  HELP_URL = 'helpUrl',
  REFERENCE = 'reference',
  LOGIN = 'author?.login',
  SAMPLE = 'sample',
  ENVIRONMENT = 'environment',
  TAGS = 'tags'
};

const obligatoryFuncProps = ['name', 'description', 'helpUrl', 'language'];

const functionPropsLabels = (key: FUNC_PROPS_FIELDS) => {
  switch (key) {
  case FUNC_PROPS_FIELDS.NAME: return 'Name';
  case FUNC_PROPS_FIELDS.DESCRIPTION: return 'Description';
  case FUNC_PROPS_FIELDS.LOGIN: return 'Author';
  case FUNC_PROPS_FIELDS.HELP_URL: return 'Help URL';
  case FUNC_PROPS_FIELDS.LANGUAGE: return 'Language';
  case FUNC_PROPS_FIELDS.REFERENCE: return 'Reference';
  case FUNC_PROPS_FIELDS.SAMPLE: return 'Sample';
  case FUNC_PROPS_FIELDS.ENVIRONMENT: return 'Environment';
  case FUNC_PROPS_FIELDS.TAGS: return 'Tags';
  }
};

const functionPropsCode = (key: string) => {
  switch (key) {
  case FUNC_PROPS_FIELDS.NAME: return 'name';
  case FUNC_PROPS_FIELDS.DESCRIPTION: return 'description';
  case FUNC_PROPS_FIELDS.LOGIN: return 'author';
  case FUNC_PROPS_FIELDS.HELP_URL: return 'helpUrl';
  case FUNC_PROPS_FIELDS.LANGUAGE: return 'language';
  case FUNC_PROPS_FIELDS.REFERENCE: return 'reference';
  case FUNC_PROPS_FIELDS.SAMPLE: return 'sample';
  case FUNC_PROPS_FIELDS.ENVIRONMENT: return 'environment';
  case FUNC_PROPS_FIELDS.TAGS: return 'tags';
  }
}; ;

const funcParamTypes = [
  DG.TYPE.BOOL,
  DG.TYPE.COLUMN_LIST,
  DG.TYPE.COLUMN,
  DG.TYPE.DATA_FRAME,
  DG.TYPE.DATE_TIME,
  DG.TYPE.FLOAT,
  DG.TYPE.GRAPHICS,
  DG.TYPE.INT,
  DG.TYPE.STRING,
];

const enum FUNC_PARAM_FIELDS {
  DIRECTION = 'direction',
  TYPE = 'propertyType',
  NAME = 'name',
  DEFAULT_VALUE = 'defaultValue',
  DESCRIPTION = 'description',
  CATEGORY = 'category',
};

const enum DIRECTION {
  INPUT= 'input',
  OUTPUT= 'output',
}

const functionParamsMapping = {
  [FUNC_PARAM_FIELDS.DIRECTION]: 'Direction',
  [FUNC_PARAM_FIELDS.TYPE]: 'Type',
  [FUNC_PARAM_FIELDS.NAME]: 'Name',
  [FUNC_PARAM_FIELDS.DEFAULT_VALUE]: 'Default value',
  [FUNC_PARAM_FIELDS.DESCRIPTION]: 'Description',
  [FUNC_PARAM_FIELDS.CATEGORY]: 'Category',
  'Direction': FUNC_PARAM_FIELDS.DIRECTION,
  'Type': FUNC_PARAM_FIELDS.TYPE,
  'Name': FUNC_PARAM_FIELDS.NAME,
  'Default value': FUNC_PARAM_FIELDS.DEFAULT_VALUE,
  'Description': FUNC_PARAM_FIELDS.DESCRIPTION,
  'Category': FUNC_PARAM_FIELDS.CATEGORY,
};

enum COMMON_TAG_NAME {
  VALIDATORS = 'validators',
  CAPTION = 'caption',
  POSTFIX = 'postfix',
  UNITS = 'units',
  EDITOR = 'editor',
}

enum OPTIONAL_TAG_NAME {
  COLUMNS = 'columns',
  TYPE = 'type',
  FORMAT = 'format',
  ALLOW_NULLS = 'allowNulls',
  ACTION = 'action',
  CHOICES = 'choices',
  SUGGESTIONS = 'suggestions',
  SEM_TYPE = 'semType',
  MIN = 'min',
  MAX = 'max',
}

type DF_TAG_NAME = COMMON_TAG_NAME | OPTIONAL_TAG_NAME.COLUMNS;
type COLUMN_TAG_NAME = COMMON_TAG_NAME |
                        OPTIONAL_TAG_NAME.COLUMNS | OPTIONAL_TAG_NAME.FORMAT |
                        OPTIONAL_TAG_NAME.ALLOW_NULLS | OPTIONAL_TAG_NAME.ACTION | OPTIONAL_TAG_NAME.SEM_TYPE

type STRING_TAG_NAME = COMMON_TAG_NAME |
                        OPTIONAL_TAG_NAME.CHOICES | OPTIONAL_TAG_NAME.SUGGESTIONS

type INT_TAG_NAME = COMMON_TAG_NAME | OPTIONAL_TAG_NAME.MIN | OPTIONAL_TAG_NAME.MAX

const COMMON_TAG_NAMES = [...Object.values(COMMON_TAG_NAME)];
const DF_TAG_NAMES = [
  ...Object.values(COMMON_TAG_NAME),
  OPTIONAL_TAG_NAME.COLUMNS,
  OPTIONAL_TAG_NAME.ACTION];
const COLUMN_TAG_NAMES = [...Object.values(COMMON_TAG_NAME), OPTIONAL_TAG_NAME.TYPE, OPTIONAL_TAG_NAME.FORMAT,
  OPTIONAL_TAG_NAME.ALLOW_NULLS, OPTIONAL_TAG_NAME.ACTION, OPTIONAL_TAG_NAME.SEM_TYPE];
const STRING_TAG_NAMES = [...Object.values(COMMON_TAG_NAME),
  /*OPTIONAL_TAG_NAME.CHOICES ,*/ OPTIONAL_TAG_NAME.SUGGESTIONS];
const INT_TAG_NAMES = [...Object.values(COMMON_TAG_NAME), OPTIONAL_TAG_NAME.MIN, OPTIONAL_TAG_NAME.MAX];

type FuncParamBase = {
  direction: DIRECTION,
  name?: string,
  defaultValue?: string,
  description?: string
}

type FuncParam =
  FuncParamBase & {
  type: DG.TYPE.BOOL | DG.TYPE.DATE_TIME | DG.TYPE.FLOAT | DG.TYPE.GRAPHICS
  options?: {tag: COMMON_TAG_NAME, value: string}[]
} | FuncParamBase & {
  type: DG.TYPE.INT
  options?: {tag: INT_TAG_NAME, value: string}[]
} | FuncParamBase & {
  type: DG.TYPE.DATA_FRAME,
  options?: {tag: DF_TAG_NAME, value: string}[]
} | FuncParamBase & {
  type: DG.TYPE.COLUMN_LIST | DG.TYPE.COLUMN,
  options?: {tag: COLUMN_TAG_NAME, value: string}[]
} | FuncParamBase & {
  type: DG.TYPE.STRING
  options?: {tag: STRING_TAG_NAME, value: string}[]
}

const optionTags = ((param: DG.Property) => {
  switch (param.propertyType) {
  case DG.TYPE.INT:
    return INT_TAG_NAMES;
  case DG.TYPE.DATA_FRAME:
    return DF_TAG_NAMES;
  case DG.TYPE.COLUMN_LIST:
  case DG.TYPE.COLUMN:
    return COLUMN_TAG_NAMES;
  case DG.TYPE.STRING:
    return STRING_TAG_NAMES;
  default:
    return COMMON_TAG_NAMES;
  }
});

enum LANGUAGE {
  JS = 'javascript', 
  PYTHON = 'python', 
  R = 'r', 
  JULIA = 'julia', 
  OCTAVE = 'octave',
  NODEJS = 'nodejs',
  GROK = 'grok'
}
const languages = ['javascript', 'python', 'r', 'julia', 'octave', 'nodejs', 'grok']
const commentSign = (lang: LANGUAGE) => {
  switch(lang) {
    case LANGUAGE.JS:
    case LANGUAGE.NODEJS:
      return '//';
    case LANGUAGE.R:
    case LANGUAGE.GROK:
    case LANGUAGE.JULIA:
    case LANGUAGE.PYTHON:
    case LANGUAGE.OCTAVE:
      return '#'
  }
}

function openFse(v: DG.View, functionCode: string) {
    const inputScriptCopy = DG.Script.create(functionCode);

    const editorView = DG.View.create();
    editorView.name = v.name;

    const openScript = () => {
      editorView.close();
      grok.shell.addView(v);
      // @ts-ignore
      let editor = v.root.lastChild.lastChild.CodeMirror;
      let doc = editor.getDoc();
      doc.setValue(myCM.getDoc().getValue());
    }
  
    const updateFuncPropValue = (propName: string, v: any) => {
      (inputScriptCopy as any)[propName] = v;
      refreshPreview();
    };

    const generateParamLine = (param: DG.Property, direction: string) => {
      const optionTagsPreview = optionTags(param)
        .map((tag) => {
          return {tag: tag, val: param.options[tag]};
        })
        .filter((value) => !!value.val)
        .map(({tag, val}) => `${tag}: ${val}`)
        .concat(...(param.category && param.category!== DEFAULT_CATEGORY)? [`category: ${param.category}`]:[])
        .join('; ');
      const result =`${commentSign(inputScriptCopy.language as LANGUAGE)}${direction}: ${param.propertyType} ${param.name? `${param.name} ` : ''}${param.defaultValue ? `= ${param.defaultValue} ` : ''}${!!optionTagsPreview.length ? `{${optionTagsPreview}} ` : ''}${param.description ? `[${param.description}]` : ''}\n`;
      return result;
    };
  
    const functionProps = {
      [functionPropsLabels(FUNC_PROPS_FIELDS.NAME)]: DG.Property.create(
        FUNC_PROPS_FIELDS.NAME, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.NAME],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.NAME, v),
        '',
      ),
      [functionPropsLabels(FUNC_PROPS_FIELDS.DESCRIPTION)]: DG.Property.create(
        FUNC_PROPS_FIELDS.DESCRIPTION, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.DESCRIPTION],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.DESCRIPTION, v),
        '',
      ),
      [functionPropsLabels(FUNC_PROPS_FIELDS.LANGUAGE)]: (() => {
        const temp = DG.Property.create(
          FUNC_PROPS_FIELDS.LANGUAGE, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.LANGUAGE],
          (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.LANGUAGE, v),
          '');
        temp.fromOptions({choices: languages});
        return temp;
      })(),
      [functionPropsLabels(FUNC_PROPS_FIELDS.HELP_URL)]: DG.Property.create(
        FUNC_PROPS_FIELDS.HELP_URL, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.HELP_URL],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.HELP_URL, v),
        '',
      ),
      [functionPropsLabels(FUNC_PROPS_FIELDS.REFERENCE)]: DG.Property.create(
        FUNC_PROPS_FIELDS.REFERENCE, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.REFERENCE],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.REFERENCE, v),
        '',
      ),
      [functionPropsLabels(FUNC_PROPS_FIELDS.LOGIN)]: DG.Property.create(
        FUNC_PROPS_FIELDS.LOGIN, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.LOGIN],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.LOGIN, v),
        '',
      ),
      [functionPropsLabels(FUNC_PROPS_FIELDS.SAMPLE)]: DG.Property.create(
        FUNC_PROPS_FIELDS.SAMPLE, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.SAMPLE],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.SAMPLE, v),
        '',
      ),
      [functionPropsLabels(FUNC_PROPS_FIELDS.ENVIRONMENT)]: DG.Property.create(
        FUNC_PROPS_FIELDS.ENVIRONMENT, DG.TYPE.STRING, (x: any) => x[FUNC_PROPS_FIELDS.ENVIRONMENT],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.ENVIRONMENT, v),
        '',
      ),
      [functionPropsLabels(FUNC_PROPS_FIELDS.TAGS)]: DG.Property.create(
        FUNC_PROPS_FIELDS.TAGS, DG.TYPE.LIST,
        (x: any) => x[FUNC_PROPS_FIELDS.TAGS],
        (x: any, v) => updateFuncPropValue(FUNC_PROPS_FIELDS.TAGS, v),
        [],
      ),
    };
  
    const getNewProps = () => {
      const newProps = [];
      for (const prop of Object.values(functionProps)) {
        if (!prop.get(inputScriptCopy) && !(inputScriptCopy as any)[prop.name]) {
          newProps.push(prop);
        }
      }
      return newProps;
    };
  
    let addButton = ui.div();
    let propsForm = ui.div();
    let inputs: HTMLElement = ui.div();
  
  
    const functionPropsInput = (prop: DG.Property) => {
      switch (prop.name) {
      case FUNC_PROPS_FIELDS.LANGUAGE:
        return ui.choiceInput(
          functionPropsLabels(prop.name as FUNC_PROPS_FIELDS),
          prop.get(inputScriptCopy) || (inputScriptCopy as any)[prop.name],
          prop.choices,
        );
      // case FUNC_PROPS_FIELDS.TAGS:
      //   return ui.multiChoiceInput(
      //     functionPropsLabels(prop.name as FUNC_PROPS_FIELDS),
      //     prop.get(inputScriptCopy) || (inputScriptCopy as any)[prop.name],
      //     prop.choices,
      //   );
      default:
        return ui.stringInput(functionPropsLabels(prop.name as FUNC_PROPS_FIELDS),
          prop.get(inputScriptCopy) || (inputScriptCopy as any)[prop.name]);
      }
    };
  
    const renderAddPropButton = () => {
      if (!getNewProps().length) return ui.div();
  
      const onItemClick = (item: DG.Property) => {
        item.set(inputScriptCopy, ' ');
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
        (prop) => menu.item(functionPropsLabels(prop.name as FUNC_PROPS_FIELDS), () => onItemClick(prop)),
      );
  
      const shevron = ui.iconFA('chevron-down');
      const button = ui.button([ui.icons.add(()=>{}), shevron], () => {
        menu.show();
      });
  
      const div = ui.div([button]);
      div.style.cssText += 'margin-left: 152px';
      return div;
    };
  
    const addFullWidthInput = (input: DG.InputBase, prop: DG.Property) => {
      (input.root.lastChild as HTMLElement).style.cssText+='width: 400px; max-width: inherit;';
      input.onInput(() => {
        prop.set(inputScriptCopy, input.stringValue);
        refreshPreview();
      });
      (input.root as HTMLInputElement).placeholder = 'Enter your value...';
      input.root.appendChild(ui.button(ui.icons.delete(()=>{}), () => {
        prop.set(inputScriptCopy, prop.defaultValue);
        refreshPreview();
        input.root.remove();
        const newAddButton = renderAddPropButton();
        addButton.replaceWith(newAddButton);
        addButton = newAddButton;
      }));
      (input.root.lastChild as HTMLElement).style.cssText+='display: inline-flex; justify-content: center; flex-direction: column';
      if (obligatoryFuncProps.includes(prop.name)) {
        input.root.style.cssText += 'padding-right: 35.375px';
        (input.root.lastChild as HTMLElement).style.cssText+='display: none;';
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
  
      const result = ui.input.form(param, optionalFuncParamsProps.filter(
        (prop) => optionTags(param).includes(prop.name as OPTIONAL_TAG_NAME)),
      );
  
      grok.shell.o = ui.divV([ui.h1(`Param: ${paramName}`), ui.block75([result])]);
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
          (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.DIRECTION, v), '');
        temp.fromOptions({choices: [DIRECTION.INPUT, DIRECTION.OUTPUT]});
        return temp;
      })(),
      DG.Property.create(FUNC_PARAM_FIELDS.NAME, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.NAME],
        (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.NAME, v), ''),
      (() => {
        const temp = DG.Property.create(FUNC_PARAM_FIELDS.TYPE, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.TYPE],
          (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.TYPE, v), '');
        temp.fromOptions({choices: funcParamTypes});
        return temp;
      })(),
      DG.Property.create(FUNC_PARAM_FIELDS.DEFAULT_VALUE, DG.TYPE.STRING,
        (x: any) => String(DG.toJs(x)[FUNC_PARAM_FIELDS.DEFAULT_VALUE] || ''),
        (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.DEFAULT_VALUE, v), ''),
      (() => {
        const temp = DG.Property.create(FUNC_PARAM_FIELDS.DESCRIPTION, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.DESCRIPTION], (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.DESCRIPTION, v), '');
        temp.fromOptions({editor: 'textarea'});
        return temp;
      })(),
      DG.Property.create(FUNC_PARAM_FIELDS.CATEGORY, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.CATEGORY],
        (x: any, v) => updateValue(x, FUNC_PARAM_FIELDS.CATEGORY, v), ''),
    ];
  
    const optionalFuncParamsProps: DG.Property[] = [
      DG.Property.create(OPTIONAL_TAG_NAME.ACTION, DG.TYPE.STRING,
        (x: any) => x.options[OPTIONAL_TAG_NAME.ACTION],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.ACTION, v), ''),
  
      DG.Property.create(OPTIONAL_TAG_NAME.ALLOW_NULLS, DG.TYPE.BOOL,
        (x: any) => x.options[OPTIONAL_TAG_NAME.ALLOW_NULLS],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.ALLOW_NULLS, v), ''),
  
      DG.Property.create(OPTIONAL_TAG_NAME.CHOICES, DG.TYPE.LIST,
        (x: any) => x.options[OPTIONAL_TAG_NAME.CHOICES],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.CHOICES, v), ''),
  
      (() => {
        const temp = DG.Property.create(OPTIONAL_TAG_NAME.COLUMNS, DG.TYPE.STRING,
          (x: any) => x.options[OPTIONAL_TAG_NAME.COLUMNS],
          (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.COLUMNS, v), '');
        temp.fromOptions({choices: ['numerical', 'categorical']});
        return temp;
      })(),
  
      DG.Property.create(OPTIONAL_TAG_NAME.FORMAT, DG.TYPE.STRING,
        (x: any) => x.options[OPTIONAL_TAG_NAME.FORMAT],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.FORMAT, v), ''),
  
      DG.Property.create(OPTIONAL_TAG_NAME.MAX, DG.TYPE.INT,
        (x: any) => x.options[OPTIONAL_TAG_NAME.MAX],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.MAX, v), ''),
  
      DG.Property.create(OPTIONAL_TAG_NAME.MIN, DG.TYPE.INT,
        (x: any) => x.options[OPTIONAL_TAG_NAME.MIN],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.MIN, v), ''),
  
      DG.Property.create(OPTIONAL_TAG_NAME.SEM_TYPE, DG.TYPE.STRING,
        (x: any) => x.options[OPTIONAL_TAG_NAME.SEM_TYPE],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.SEM_TYPE, v), ''),
  
      DG.Property.create(OPTIONAL_TAG_NAME.SUGGESTIONS, DG.TYPE.STRING,
        (x: any) => x.options[OPTIONAL_TAG_NAME.SUGGESTIONS],
        (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.SUGGESTIONS, v), ''),
  
      (() => {
        const temp = DG.Property.create(OPTIONAL_TAG_NAME.TYPE, DG.TYPE.STRING,
          (x: any) => x.options[OPTIONAL_TAG_NAME.TYPE],
          (x: any, v) => updateValue(x, OPTIONAL_TAG_NAME.TYPE, v), '');
        temp.fromOptions({choices: ['numerical', 'categorical', 'dateTime']});
        return temp;
      })(),
    ];
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
    paramsDF.onCurrentRowChanged.subscribe(() => onFunctionParamClick((paramsDF.currentRow as any)['Name']));
  
    const paramsGrid = DG.Grid.create(paramsDF);
    paramsGrid.root.style.width = '100%';
    paramsGrid.dataFrame?.getCol(functionParamsMapping[FUNC_PARAM_FIELDS.TYPE as keyof typeof functionParamsMapping])
      .setTag(DG.TAGS.CHOICES, `["${funcParamTypes.join(`", "`)}"]`);
    paramsGrid.dataFrame?.getCol(functionParamsMapping[FUNC_PARAM_FIELDS.DIRECTION as keyof typeof functionParamsMapping])
      .setTag(DG.TAGS.CHOICES, `["${[DIRECTION.INPUT, DIRECTION.OUTPUT].join(`", "`)}"]`);
    paramsGrid.onCellPrepare((gc) => {
      const deleteBtn = (name: string) => ui.div(
        ui.icons.delete(() => {
          console.log('name to delete:', name);
          functionParamsCopy = functionParamsCopy.filter((param) => param.name !== name);
          functionParamsState.next(functionParamsCopy);
        }, 'Remove the param'), {style: {'text-align': 'center', 'margin': '6px'}},
      );
  
      if (gc.gridColumn.name === '+') {
        gc.gridColumn.cellType = 'html';
        if (gc.isTableCell) {
          gc.style.element =
          ui.divH([
            deleteBtn(gc.grid.dataFrame?.get('Name', gc.gridRow)),
          ]);
        }
  
        if (gc.isColHeader) {
          gc.customText = '';
        }
      }
    });
  
    const area = ui.textInput('', '');
    const myCM = CodeMirror.fromTextArea((area.input as HTMLTextAreaElement));
  
    const previewDiv = ui.panel([ui.divV([ui.h1('Code preview'), area])]);
    previewDiv.style.flexGrow = '1';
  
    propsForm = functionPropsForm();
  
    const addParamBtn = () => ui.button(
      [
        'Add parameter',
        ui.div(ui.icons.add(() => {
        }), {style: {'text-align': 'center', 'margin': '6px'}}),
      ],
      () => {
        const newParam = {
          [FUNC_PARAM_FIELDS.DIRECTION]: DIRECTION.INPUT,
          [FUNC_PARAM_FIELDS.NAME]: 'newParam',
          [FUNC_PARAM_FIELDS.TYPE]: DG.TYPE.BOOL,
          [FUNC_PARAM_FIELDS.DEFAULT_VALUE]: false,
          [FUNC_PARAM_FIELDS.DESCRIPTION]: '',
          [FUNC_PARAM_FIELDS.CATEGORY]: '',
        };
  
        const t = DG.Property.create(newParam.name, newParam.propertyType, (x: any) => x, (x: any, v) => x = v);
        t.options.direction = newParam.direction;
        t.description = newParam.description;
        functionParamsCopy.push(t);
        functionParamsState.next(functionParamsCopy);
        refreshPreview();
      },
    );
  
    const tabs = ui.tabControl({
      'PROPERTIES': propsForm,
      'PARAMETERS': ui.divV([
        addParamBtn(),
        paramsGrid.root,
      ]),
    });
    tabs.root.style.width = '100%';
    tabs.root.style.flexGrow = '3';
  
    editorView.append(
      ui.divV([
        tabs,
        previewDiv,
      ]),
    );
    editorView.box = true;
    editorView.setRibbonPanels([
      [
        ui.iconFA('eye', () => {
          previewDiv.hidden ? previewDiv.hidden = false : previewDiv.hidden = true;
        }),
        ui.iconFA('code', () => openScript(), 'Open function editor'),
      ],
    ]);
  
    const refreshPreview = () => {
      let result = '';
      Object.values(functionProps).map((propField) => {
        const propValue = propField.get(inputScriptCopy) || (inputScriptCopy as any)[propField.name];
        if (!!propValue && !!propValue.length) {
          result += `${commentSign(inputScriptCopy.language as LANGUAGE)}${functionPropsCode(propField.name as FUNC_PROPS_FIELDS)}: ${propValue}\n`;
        }
      });
      functionParamsCopy.map((param) => {
        result += generateParamLine(param, param.options.direction);
      });
      result += inputScriptCopy.script.substring(inputScriptCopy.script.indexOf('\n',inputScriptCopy.script.lastIndexOf(commentSign(inputScriptCopy.language as LANGUAGE)))+1);
      myCM.setValue(result);
    };
  
    v.close()
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
        const newParam = functionParamsCopy[functionParamsCopy.length-1];
        paramsDF.rows.addNew(
          [
            newParam.options.direction, newParam.name,
            newParam.propertyType, newParam.defaultValue,
            newParam.description, newParam.category,
          ],
        );
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
  
    paramsGrid.onCellValueEdited.subscribe((editedCell)=> {
      const rowIndex = functionParamsCopy.findIndex((param) => param.name === (editedCell.tableRow as any)['Name']);
      if (rowIndex) {
        (functionParamsCopy[rowIndex] as any)
          [functionParamsMapping[editedCell.cell.column.name as keyof typeof functionParamsMapping]] = editedCell.cell.value;
        functionParamsCopy[rowIndex].options
          [functionParamsMapping[editedCell.cell.column.name as keyof typeof functionParamsMapping]] = editedCell.cell.value;
        functionParamsState.next(functionParamsCopy);
      }
    });
}
