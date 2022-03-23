import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { caption, enabled, checkHTMLElement } from './utils';

category('UI: Forms from props', () => {
  let v: DG.View;
  const properties = [
    {
      "name": "reactorType",
      "type": "string",
      "choices": ["Experimental", "Production"],
      "caption": "My Custom Caption"
    },
    {
      "name": "structure",
      "type": "string",
      "semType": "Molecule",
      "units": '#',
    },
    {
      "name": "started",
      "type": DG.TYPE.DATE_TIME,
    },
    {
      "name": "rating",
      "type": "int",
      "editor": "slider",
      "min": 0,
      "max": 10,
      "units": 'kg',
    }
  ].map((p) => DG.Property.fromOptions(p)); 
  
  const source = {
    reactorType: 'Production',
    structure: 'CC(=O)OC1=CC=CC=C1C(=O)O',
    rating: 7,
    department: 'IT'
  };

  before(async () => {
    v = grok.shell.newView('');
    v.append(ui.input.form(source, properties))
  });

  test('input types', async () => {
    const allInputTypesExist = 
      !!v.root.querySelector('.ui-input-slider.ui-input-root') && 
      !!v.root.querySelector('.ui-input-date.ui-input-root') && 
      !!v.root.querySelector('.ui-input-text.ui-input-root') && 
      !!v.root.querySelector('.ui-input-choice.ui-input-root')

    expect(allInputTypesExist, true);
  })

  test('input captions', async () => {
    const allSpans = Array.from(v.root.querySelectorAll('span'))

    const allCaptionsExist = 
      !!allSpans.find((el) => el.innerText === 'My Custom Caption') &&
      !!allSpans.find((el) => el.innerText === 'Structure') &&
      !!allSpans.find((el) => el.innerText === 'Started') &&
      !!allSpans.find((el) => el.innerText === 'Rating');

    expect(allCaptionsExist, true);
  })

  test('input units', async () => {
    const allUnits = Array.from(v.root.querySelectorAll('label'))

    const allUnitsExist = 
      !!allUnits.find((el) => el.innerText === '#') &&
      !!allUnits.find((el) => el.innerText === 'kg');

    expect(allUnitsExist, true);
  })

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});

