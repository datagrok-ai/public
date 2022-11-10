import * as DG from 'datagrok-api/dg';
import { RGROUP_CAP_GROUP_NAME, RGROUP_CAP_GROUP_SMILES, jsonSdfMonomerLibDict, MONOMER_SYMBOL, RGROUP_ALTER_ID, RGROUPS, RGROUP_LABEL, SDF_MONOMER_NAME } from "./constants";

export function getParts(subParts: string[], s: string): string[] {
  let j = 0;
  let allParts: string[] = [];
  for (let k = 0; k < subParts.length; ++k) {
    let indexOfMonomer = s.indexOf(subParts[k]);
    let helmBeforeMonomer = s.slice(j, indexOfMonomer);
    allParts.push(helmBeforeMonomer);
    allParts.push(subParts[k]);
    s = s.substring(indexOfMonomer + subParts[k].length);
  }
  allParts.push(s);
  return allParts;
}

export function parseHelm(s: string) {
  var sections = split(s, '$');
  s = sections[0];
  var monomers = [];
  //@ts-ignore
  if (!scil.Utils.isNullOrEmpty(s)) {
    var seqs = split(s, '|');
    for (var i = 0; i < seqs.length; ++i) {
      var e = detachAnnotation(seqs[i]);
      s = e.str;

      var p = s.indexOf('{');

      s = s.substring(p + 1);
      p = s.indexOf('}');
      s = s.substring(0, p);

      var ss = split(s, '.');
      for (var monomer of ss) {
        if (monomer.includes('(') && monomer.includes(')')) {
          var elements = monomer.replace(/[()]/g, '').split('');
          for (var el of elements) {
            monomers.push(el);
          }
        } else if (monomer.includes('[') && monomer.includes(']')) {
          var element = monomer.match(/(?<=\[).*?(?=\])/g, '');
          monomers.push(element[0]);
        } else {
          monomers.push(monomer);
        }
      }
    }
  }
  return monomers;
}

export function findMonomers(helmString: string) {
  //@ts-ignore
  const types = Object.keys(org.helm.webeditor.monomerTypeList());
  const monomers: any = [];
  const monomer_names: any = [];
  for (var i = 0; i < types.length; i++) {
    //@ts-ignore
    monomers.push(new scil.helm.Monomers.getMonomerSet(types[i]));
    Object.keys(monomers[i]).forEach(k => {
      monomer_names.push(monomers[i][k].id);
    });
  }
  const split_string = parseHelm(helmString);
  return new Set(split_string.filter(val => !monomer_names.includes(val)));
}

function split(s: string, sep: string) {
  var ret = [];
  var frag = '';
  var parentheses = 0;
  var bracket = 0;
  var braces = 0;
  var quote = 0;
  for (var i = 0; i < s.length; ++i) {
    var c = s.substring(i, i + 1);
    if (c == sep && bracket == 0 && parentheses == 0 && braces == 0 && quote == 0) {
      ret.push(frag);
      frag = '';
    } else {
      frag += c;
      if (quote > 0) {
        if (c == '\\' && i + 1 < s.length) {
          ++i;
          var c2 = s.substring(i, i + 1);
          frag += c2;
          c += c2;
        }
      }
      if (c == '\"') {
        if (!(i > 0 && s.substring(i - 1, i) == '\\'))
          quote = quote == 0 ? 1 : 0;
      } else if (c == '[')
        ++bracket;
      else if (c == ']')
        --bracket;
      else if (c == '(')
        ++parentheses;
      else if (c == ')')
        --parentheses;
      else if (c == '{')
        ++braces;
      else if (c == '}')
        --braces;
    }
  }
  ret.push(frag);
  return ret;
}

function detachAnnotation(s: string) {
  var ret = _detachAppendix(s, '\"');
  if (ret.tag != null)
    return ret;

  var r = _detachAppendix(s, '\'');
  return {tag: ret.tag, repeat: r.tag, str: r.str};
}

function _detachAppendix(s: string, c: string) {
  var tag = null;
  //@ts-ignore
  if (scil.Utils.endswith(s, c)) {
    var p = s.length - 1;
    while (p > 0) {
      p = s.lastIndexOf(c, p - 1);
      if (p <= 0 || s.substring(p - 1, p) != '\\')
        break;
    }

    if (p > 0 && p < s.length - 1) {
      tag = s.substring(p + 1, s.length - 1);
      s = s.substring(0, p);
    }
  }
  if (tag != null)
    tag = tag.replace(new RegExp('\\' + c, 'g'), c);
  return {tag: unescape(tag), str: s};
}

function unescape(s: string) {
  //@ts-ignore
  if (scil.Utils.isNullOrEmpty(s))
    return s;

  return s.replace(/[\\]./g, function(m) {
    switch (m) {
    case '\\r':
      return '\r';
    case '\\n':
      return '\n';
    case '\\t':
      return '\t';
    default:
      return m.substring(1);
    }
  });
}
