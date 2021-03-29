/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

interface Rule {
  column: string;
  type: string;
  pattern: string;
  matcher?: any;
}

interface RuleSet {
  name: string;
  rules: Rule[];
}

// @ts-ignore
let ruleSet: RuleSet = {
  name: 'SDTM DM rules',
  rules: [
    {
      column: 'age',
      type: 'int',
      pattern: '18-90'
    },
    {
      column: 'height',
      type: 'numerical',
      pattern: '18-90'
    },
    {
      column: 'sex',
      type: 'string',
      pattern: 'in (F, M)'
    },
    {
      column: 'started',
      type: 'datetime',
      pattern: '> 1/1/1900'
    }
  ]
};

/** Returns error message, or null if valid. */
function validateCell(cell: any, ruleset: RuleSet) : string | null {
  let name = cell.column.name.toLowerCase();

  for (let rule of ruleset.rules)
    if (rule.column == name) {
      rule.matcher = rule.matcher ?? DG.ValueMatcher.fromColumn(cell.column, rule.pattern);
      let error = rule.matcher.validate(cell.value);
      if (error !== null)
        return error;
    }

  return null;
}

function main() : void {
  console.log('foo');
}

let cell = {column: { name: 'age' }, value: 999 };
validateCell(cell, ruleSet);

console.log('ddd');