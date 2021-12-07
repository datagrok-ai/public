import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import CodeMirror from 'codemirror';
import 'codemirror/lib/codemirror.css';
import {_package} from '../package';

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

const functionPropsCode = (key: FUNC_PROPS_FIELDS) => {
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
};

type FuncPropsState = {
  [FUNC_PROPS_FIELDS.NAME]: string | null,
  [FUNC_PROPS_FIELDS.DESCRIPTION]: string | null,
  [FUNC_PROPS_FIELDS.LANGUAGE]: string | null,
  [FUNC_PROPS_FIELDS.HELP_URL]: string | null,
  [FUNC_PROPS_FIELDS.REFERENCE]: string | null,
  [FUNC_PROPS_FIELDS.LOGIN]: string | null,
  [FUNC_PROPS_FIELDS.SAMPLE]: string | null,
  [FUNC_PROPS_FIELDS.ENVIRONMENT]: string | null,
  [FUNC_PROPS_FIELDS.TAGS]: string | null,
}

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

const functionParamsLabels = (key: FUNC_PARAM_FIELDS) => {
  switch (key) {
  case FUNC_PARAM_FIELDS.DIRECTION: return 'Direction';
  case FUNC_PARAM_FIELDS.TYPE: return 'Type';
  case FUNC_PARAM_FIELDS.NAME: return 'Name';
  case FUNC_PARAM_FIELDS.DEFAULT_VALUE: return 'Default value';
  case FUNC_PARAM_FIELDS.DESCRIPTION: return 'Description';
  case FUNC_PARAM_FIELDS.CATEGORY: return 'Category';
  }
};

enum COMMON_TAG_NAME {
  VALIDATORS = 'validators',
  CAPTION = 'caption',
  POSTFIX = 'postfix',
  UNITS = 'units'
}

enum OPTIONAL_TAG_NAME {
  COLUMNS = 'columns',
  CATEGORICAL = 'categorical',
  TYPE = 'type',
  FORMAT = 'format',
  ALLOW_NULLS = 'allowNulls',
  ACTION = 'action',
  CHOICES = 'choices',
  SUGGESTIONS = 'suggestions',
  SEM_TYPE = 'semType'
}

type DF_TAG_NAME = COMMON_TAG_NAME | OPTIONAL_TAG_NAME.COLUMNS | OPTIONAL_TAG_NAME.CATEGORICAL;
type COLUMN_TAG_NAME = COMMON_TAG_NAME |
                        OPTIONAL_TAG_NAME.TYPE | OPTIONAL_TAG_NAME.FORMAT |
                        OPTIONAL_TAG_NAME.ALLOW_NULLS | OPTIONAL_TAG_NAME.ACTION | OPTIONAL_TAG_NAME.SEM_TYPE

type STRING_TAG_NAME = COMMON_TAG_NAME |
                        OPTIONAL_TAG_NAME.CHOICES | OPTIONAL_TAG_NAME.SUGGESTIONS

const DF_TAG_NAMES = [...Object.values(COMMON_TAG_NAME), OPTIONAL_TAG_NAME.COLUMNS, OPTIONAL_TAG_NAME.CATEGORICAL];
const COLUMN_TAG_NAMES = [...Object.values(COMMON_TAG_NAME), OPTIONAL_TAG_NAME.TYPE, OPTIONAL_TAG_NAME.FORMAT,
  OPTIONAL_TAG_NAME.ALLOW_NULLS, OPTIONAL_TAG_NAME.ACTION, OPTIONAL_TAG_NAME.SEM_TYPE];
const STRING_TAG_NAMES = [...Object.values(COMMON_TAG_NAME), OPTIONAL_TAG_NAME.CHOICES, OPTIONAL_TAG_NAME.SUGGESTIONS];

type FuncParamBase = {
  direction: DIRECTION,
  name?: string,
  defaultValue?: string,
  description?: string
}

type FuncParam =
    FuncParamBase & {
  type: DG.TYPE.BOOL | DG.TYPE.DATE_TIME | DG.TYPE.FLOAT | DG.TYPE.GRAPHICS | DG.TYPE.INT
  options?: {tag: COMMON_TAG_NAME, value: string}[]
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

const generateParamLine = (param: DG.Property, direction: string) => {
  const optionTags = `{${
    COLUMN_TAG_NAMES
      .map((tag) => {
        console.log(param.name, tag, param.options(), (param.options() as any)[tag]);
        return {tag: tag, val: (param.options() as any)[tag]};
      })
      .filter((value) => !!value.val)
      .map(({tag, val}) => `${tag}:${val}`)
      .join('; ')
  }}`;
  const result =`#${direction}: ${param.propertyType} ${param.name || ''} ${param.defaultValue ? `= ${param.defaultValue}` : ''} ${optionTags} ${param.description ? ` [${param.description}]` : ''}\n`;
  return result;
};

export function _functionEditor(functionCode: string) {
  const inputScript = DG.Script.create(functionCode);

  const editorView = DG.View.create();
  editorView.name = 'Function Editor';

  const functionPropsCopy: FuncPropsState = {
    [FUNC_PROPS_FIELDS.NAME]: inputScript.name || null,
    [FUNC_PROPS_FIELDS.DESCRIPTION]: inputScript.description || null,
    [FUNC_PROPS_FIELDS.LANGUAGE]: inputScript.language || null,
    [FUNC_PROPS_FIELDS.HELP_URL]: inputScript.helpUrl || null,
    [FUNC_PROPS_FIELDS.REFERENCE]: inputScript.reference || null,
    [FUNC_PROPS_FIELDS.LOGIN]: inputScript.author?.login || null,
    [FUNC_PROPS_FIELDS.SAMPLE]: inputScript.sample || null,
    [FUNC_PROPS_FIELDS.ENVIRONMENT]: inputScript.environment || null,
    [FUNC_PROPS_FIELDS.TAGS]: inputScript.tags?.join(' ') || null,
  };

  const getNewProps = () => {
    const newProps = [];
    for (const key of Object.keys(functionPropsCopy)) {
      if (functionPropsCopy[key as keyof FuncPropsState] == null) {
        newProps.push(key);
      }
    }
    return newProps;
  };

  let addButton = ui.div();
  let propsForm = ui.div();
  let inputs = ui.div();

  const renderAddPropButton = () => {
    if (!getNewProps().length) return ui.div();

    const onItemClick = (item: string) => {
      functionPropsCopy[item as keyof FuncPropsState] = '';
      const newInputs = ui.inputs(
        Object.keys(functionPropsCopy)
          .filter((prop) => functionPropsCopy[prop as keyof FuncPropsState] !== null)
          .map((prop) => addFullWidthInput(
            ui.stringInput(
              functionPropsLabels(prop as keyof FuncPropsState),
              functionPropsCopy[prop as keyof FuncPropsState] || ''), prop,
          )),
      );
      inputs.replaceWith(newInputs);
      inputs = newInputs;
      refreshPreview();
      const newAddButton = renderAddPropButton();
      addButton.replaceWith(newAddButton);
      addButton = newAddButton;
    };
    const menu = DG.Menu.popup().items(getNewProps(), onItemClick);

    const shevron = ui.iconFA('chevron-down');
    const button = ui.button([ui.icons.add(()=>{}), shevron], () => {
      menu.show();
    });

    const div = ui.div([button]);
    div.style.cssText += 'margin-left: 152px';
    return div;
  };

  const addFullWidthInput = (input: DG.InputBase, propName: string) => {
    (input.root.lastChild as HTMLElement).style.cssText+='width: 400px; max-width: inherit;';
    input.onInput(() => {
      functionPropsCopy[propName as keyof typeof functionPropsCopy] = input.stringValue;
      refreshPreview();
    });
    (input.root as HTMLInputElement).placeholder = 'Enter your value...';
    input.root.append(ui.button(ui.icons.delete(()=>{}), () => {
      functionPropsCopy[propName as keyof typeof functionPropsCopy] = null;
      refreshPreview();
      input.root.remove();
      const newAddButton = renderAddPropButton();
      addButton.replaceWith(newAddButton);
      addButton = newAddButton;
    }));
    (input.root.lastChild as HTMLElement).style.cssText+='display: inline-flex; justify-content: center; flex-direction: column';
    if (obligatoryFuncProps.includes(propName)) {
      input.root.style.cssText += 'padding-right: 35.375px';
      (input.root.lastChild as HTMLElement).style.cssText+='display: none;';
    }
    return input;
  };

  const functionPropsForm = () => {
    addButton = renderAddPropButton();

    inputs = ui.inputs(
      Object.keys(functionPropsCopy)
        .filter((prop) => !!functionPropsCopy[prop as keyof FuncPropsState])
        .map((prop) => addFullWidthInput(
          ui.stringInput(
            functionPropsLabels(prop as keyof FuncPropsState),
            functionPropsCopy[prop as keyof FuncPropsState] || ''), prop,
        )),
    );
    return ui.panel([
      ui.divV([
        inputs,
        addButton,
      ]),
    ]);
  };

  let functionInputParamsCopy = inputScript.inputs;
  let functionOutputParamsCopy = inputScript.outputs;

  const iconsByTypes = (type: DG.TYPE) => {
    switch (type) {
    case DG.TYPE.BOOL: return ui.iconImage(DG.TYPE.BOOL, `${_package.webRoot}/icons/Bool.svg`);
    case DG.TYPE.COLUMN_LIST: return ui.iconImage(DG.TYPE.BOOL, `${_package.webRoot}/icons/Column-list.svg`);
    case DG.TYPE.COLUMN: return ui.iconImage(DG.TYPE.COLUMN, `${_package.webRoot}/icons/Column.svg`);
    case DG.TYPE.DATA_FRAME: return ui.iconImage(DG.TYPE.DATA_FRAME, `${_package.webRoot}/icons/Dataframe.svg`);
    case DG.TYPE.DATE_TIME: return ui.iconImage(DG.TYPE.COLUMN, `${_package.webRoot}/icons/Datetime.svg`);
    case DG.TYPE.FLOAT: return ui.iconImage(DG.TYPE.FLOAT, `${_package.webRoot}/icons/Double.svg`);
    case DG.TYPE.GRAPHICS: return ui.iconImage(DG.TYPE.GRAPHICS, `${_package.webRoot}/icons/Graphics.svg`);
    case DG.TYPE.INT: return ui.iconImage(DG.TYPE.FLOAT, `${_package.webRoot}/icons/Int.svg`);
    case DG.TYPE.STRING: return ui.iconImage(DG.TYPE.STRING, `${_package.webRoot}/icons/String.svg`);
    default: return ui.iconImage('Blob', `${_package.webRoot}/icons/Blob.svg`);
    }
  };

  const onFunctionParamClick = (paramName: string, card: HTMLElement, direction: DIRECTION) => {
    const acc = ui.accordion();
    acc.root.style.cssText += 'min-width: 400px';

    acc.addPane('General', ()=> {
      let paramIndex: number;
      let param: DG.Property;

      if (direction === DIRECTION.INPUT) {
        paramIndex = functionInputParamsCopy.findIndex((input) => input.name === paramName);
        param = functionInputParamsCopy[paramIndex];
      } else {
        paramIndex= functionOutputParamsCopy.findIndex((input) => input.name === paramName);
        param = functionOutputParamsCopy[paramIndex];
      }

      if (!param) return ui.div('');

      const updateCard = (propName: string, v: any) => {
        (param as any)[propName] = v;
        const newCard = renderCard(param, direction);
        card.replaceWith(newCard);
        card = newCard;
        refreshPreview();
      };

      const funcParams: DG.Property[] = [
        (() => {
          const temp = DG.Property.create(FUNC_PARAM_FIELDS.DIRECTION, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.DIRECTION], (x: any, v) => updateCard(FUNC_PARAM_FIELDS.DIRECTION, v), '');
          temp.options({choices: ['input', 'output']});
          return temp;
        })(),
        (() => {
          const temp = DG.Property.create(FUNC_PARAM_FIELDS.TYPE, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.TYPE], (x: any, v) => updateCard(FUNC_PARAM_FIELDS.TYPE, v), '');
          temp.options({choices: funcParamTypes});
          return temp;
        })(),
        DG.Property.create(FUNC_PARAM_FIELDS.NAME, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.NAME], (x: any, v) => updateCard(FUNC_PARAM_FIELDS.NAME, v), ''),
        DG.Property.create(FUNC_PARAM_FIELDS.DEFAULT_VALUE, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.DEFAULT_VALUE], (x: any, v) => updateCard(FUNC_PARAM_FIELDS.DEFAULT_VALUE, v), ''),
        (() => {
          const temp = DG.Property.create(FUNC_PARAM_FIELDS.DESCRIPTION, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.DESCRIPTION], (x: any, v) => updateCard(FUNC_PARAM_FIELDS.DESCRIPTION, v), '');
          temp.options({editor: 'textarea'});
          return temp;
        })(),
        DG.Property.create(FUNC_PARAM_FIELDS.CATEGORY, DG.TYPE.STRING, (x: any) => x[FUNC_PARAM_FIELDS.CATEGORY], (x: any, v) => updateCard(FUNC_PARAM_FIELDS.CATEGORY, v), ''),
      ];

      const paramProps = funcParams
        .filter((param) => param.name !== FUNC_PARAM_FIELDS.DIRECTION);

      const result = ui.input.form(param, paramProps);

      return ui.div(result);
    });

    acc.addPane('Tags', ()=> {
      return ui.div();
    });

    grok.shell.o = ui.divV([ui.h1(`Param: ${paramName}`), acc.root]);
  };

  const funcParams: DG.Property[] = [
    DG.Property.jsString(FUNC_PARAM_FIELDS.DIRECTION),
    DG.Property.jsString(FUNC_PARAM_FIELDS.TYPE),
    DG.Property.jsString(FUNC_PARAM_FIELDS.NAME),
    DG.Property.jsString(FUNC_PARAM_FIELDS.DEFAULT_VALUE),
    DG.Property.jsString(FUNC_PARAM_FIELDS.DESCRIPTION),
    DG.Property.jsString(FUNC_PARAM_FIELDS.CATEGORY),
  ];


  const renderCard = (functionParam: DG.Property, direction: DIRECTION) => {
    const iconStyle = {style: {
      'display': 'flex',
      'flex-direction': 'column',
      'justify-content': 'center',
      'width': '32.5px',
      'align-items': 'center',
    }};

    const editButtonDiv = ui.div(ui.button(ui.icons.edit(()=>{}), () =>{
      onFunctionParamClick(functionParam.name, tempCard, direction);
    }), {
      style: {
        ...iconStyle.style,
        'visibility': 'hidden',
      },
    });

    const tempCard = ui.divH([
      ui.div(iconsByTypes(functionParam.propertyType), iconStyle),
      ui.divV([
        ui.span([functionParam.name], {style: {'font-weight': 'bold'}}),
        ui.span([funcParams
          .filter((param) => param.name !== FUNC_PARAM_FIELDS.NAME)
          .map((param) => (functionParam as any)[param.name])
          .filter((paramValue) => !!paramValue).join(' | '),
        ]),
      ], {style: {'width': '100%', 'padding': '5px'}}),
      editButtonDiv,
      ui.div(ui.button(ui.icons.delete(()=>{}), () =>{
        if (direction === DIRECTION.INPUT) {
          functionInputParamsCopy = functionInputParamsCopy.filter((param) => param !== functionParam);
        } else {
          functionOutputParamsCopy = functionOutputParamsCopy.filter((param) => param !== functionParam);
        }
        refreshPreview();
        tempCard.remove();
      }), iconStyle),
    ], {style: {'border-bottom': 'thin solid #ECEFF2'}});
    tempCard.addEventListener('click', () => onFunctionParamClick(functionParam.name, tempCard, direction));
    tempCard.addEventListener('mouseover', () => editButtonDiv.style.visibility = 'visible');
    tempCard.addEventListener('mouseout', () => editButtonDiv.style.visibility = 'hidden');
    return tempCard;
  };

  let inputParamCards = ui.divV(
    functionInputParamsCopy.map((functionParam) => renderCard(functionParam, DIRECTION.INPUT)),
  );
  const outputParamCards = ui.divV(
    functionOutputParamsCopy.map((functionParam) => renderCard(functionParam, DIRECTION.OUTPUT)),
  );

  const area = ui.textInput('', '');
  const myCM = CodeMirror.fromTextArea((area.input as HTMLTextAreaElement));

  const previewDiv = ui.divV([ui.h1('Code preview'), area]);
  previewDiv.style.flexGrow = '1';

  propsForm = functionPropsForm();

  const renderAddParamButton = () => {
    const paramProps = funcParams
      .filter((param) => param.name !== FUNC_PARAM_FIELDS.DIRECTION);

    const newParam: DG.PropertyOptions = {};

    const onItemClick = () => {
      ui.dialog('Add parameter')
        .add(ui.input.form(newParam, paramProps))
        .onOK(() => {
          functionInputParamsCopy.push(DG.Property.fromOptions(newParam));
          const newInputParamsCards = ui.divV(
            functionInputParamsCopy.map((functionParam) => renderCard(functionParam, DIRECTION.INPUT)),
          );
          inputParamCards.replaceWith(newInputParamsCards);
          inputParamCards = newInputParamsCards;
          refreshPreview();
        })
        .show();
    };

    const menu = DG.Menu.popup().items(funcParamTypes, onItemClick);

    const shevron = ui.iconFA('chevron-down');
    const button = ui.button([ui.icons.add(()=>{}), shevron], () => {
      menu.show();
    });

    const div = ui.div([button]);
    return div;
  };
  const addInputButton = renderAddParamButton();
  const addOutputButton = renderAddParamButton();

  const tabs = ui.tabControl({
    'PROPERTIES': propsForm,
    'PARAMETERS': ui.divV([
      ui.div(ui.h1('Input: '), {style: {'padding': '12px 12px 0px 12px'}}),
      inputParamCards,
      addInputButton,
      ui.div(ui.h1('Output: '), {style: {'padding': '12px 12px 0px 12px'}}),
      outputParamCards,
      addOutputButton,
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
      ui.bigButton('Save', () => grok.shell.info('plus')),
    ],
  ]);

  const refreshPreview = () => {
    let result = '';
    Object.keys(functionPropsCopy).map((propField) => {
      const propValue = functionPropsCopy[propField as keyof FuncPropsState];
      if (propValue) result += `#${functionPropsCode(propField as keyof FuncPropsState)}: ${propValue}\n`;
    });
    functionInputParamsCopy.map((param) => {
      result += generateParamLine(param, DIRECTION.INPUT);
    });
    functionOutputParamsCopy.map((param) => {
      result += generateParamLine(param, DIRECTION.OUTPUT);
    });
    myCM.setValue(result);
  };

  grok.shell.addView(editorView);
  refreshPreview();
}
