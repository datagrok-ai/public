---
title: "Text functions"
---

As parameters of the function, you can pass a regular string or a column name. To pass a regular string, it must be
enclosed in double quotes, like `"my string"` or in apostrophes, like  `'my string'`. To pass a column cell, you can use
the syntax `${columnName}`.

> Note that character and substring indexing in strings and lists starts at 0.

*Function List:*

- [Add](#add)
- [Contains](#contains)
- [EndsWith](#endswith)
- [Eq](#eq)
- [IsEmpty](#isempty)
- [IsNotEmpty](#isnotempty)
- [Length](#length)
- [NotEq](#noteq)
- [ParseFloat](#parsefloat)
- [ParseInt](#parseint)
- [RegExpContains](#regexpcontains)
- [RegExpExtract](#regexpextract)
- [RegExpReplace](#rexpreplace)
- [ReplaceAll](#replaceall)
- [SplitString](#splitstring)
- [StartsWith](#startswith)
- [StrFind](#strfind)
- [StrLeft](#strleft)
- [StrRight](#strright)
- [StrRepeat](#strrepeat)
- [Substring](#substring)
- [ToLowerCase](#tolowercase)
- [ToUpperCase](#touppercase)
- [Trim](#trim)

## <a name="add"></a>Add(`s1`, `s2`)

Append the string `s2` at the end of the string `s1` and returns the result obtained.

```javascript
Add("Police", "man")    // "Policeman"
```

## <a name="contains"></a>Contains(`s`, `sub`)

Checks if the string `s` contains a match of `sub`.

```javascript
Contains("Stormy Weather", "Sunny")      // false
Contains("Stormy Weather", "Weather")    // true
```

## <a name="endswith"></a>EndsWith(`s`, `postfix`)

Checks if the string `s` ends with a match of `postfix`.

```javascript
EndsWith("White Christmas", "White")        // false
EndsWith("White Christmas", "Christmas")    // true
```

## <a name="eq"></a>Eq(`s1`, `s2`)

Returns true if the string `s1` equal to `s2` and false otherwise.

```javascript
Eq("Sky", "Sky")    // true
```

## <a name="isempty"></a>IsEmpty(`s`)

Returns true if the string `s` is empty or null, and false otherwise.

```javascript
IsEmpty("")         // true
IsEmpty(null)       // true
IsEmpty("Green")    // false
```

## <a name="isnotempty"></a>IsNotEmpty(`s`)

Returns true if the string `s` is not empty and not null, and false otherwise.

```javascript
IsNotEmpty("")         // false
IsNotEmpty(null)       // false
IsNotEmpty("Green")    // true
```

## <a name="length"></a>Length(`s`)

Returns the length of the string `s`.

```javascript
Length("Text")    // 4
```

## <a name="noteq"></a>NotEq(`x`, `y`)

Returns false if the string `s1` equal to `s2` and true otherwise.

```javascript
Eq("Sky", "Sky")    // false
```

## <a name="parsefloat"></a>ParseFloat(`s`)

Parse `s` as a, possibly signed, real literal and return its value.

```javascript
ParseFloat("2025")      // -2025
ParseFloat("-012.78")   // -12.78
```

## <a name="parseint"></a>ParseInt(`s`)

Parse `s` as a, possibly signed, integer literal and return its value.

```javascript
ParseInt("2025")    // 2025
ParseInt("-012")    // -12
```

## <a name="regexpcontains"></a>RegExpContains(`s`, `pattern`)

Checks if the string `s` contains a string matches a regular expression `pattern`.

```javascript
RegExpContains("Stormy Weather", "Sunny")      // false
RegExpContains("name@gmail.com", "(\W|^)[\w.\-]{0,25}@(hotmail|gmail)\.com(\W|$)")    // true
```

## <a name="regexpextract"></a>RegExpExtract(`s`, `pattern`, `n`)

Returns the `n`-th part of a string `s` that matches a regular expression `pattern`.

```javascript
RegExpExtract("Hello World!", "l+", 0)    // "ll"
RegExpExtract("Hello World!", "l+", 1)    // "l"
```

## <a name="rexpreplace"></a>RegExpReplace(`s`, `pattern`, `sub`)

Returns the string after replacing a substring `sub` in string `s` according to a regular expression `pattern`.

```javascript
RegExpReplace("Hello World!", "l+", "LL")    // "HeLLo WorLLd!"
```

## <a name="replaceall"></a>ReplaceAll(`s`, `from`, `replace`)

Replaces all substrings in string `s` that match `from` with `replace` and returns the result obtained.

```javascript
ReplaceAll("New York", "York", "Orleans")    // "New Orleans"
ReplaceAll("Every", "", ".")                 // ".E.v.e.r.y."
ReplaceAll("moto", "o", "")                  // "mt"
```

## <a name="splitstring"></a>SplitString(`s`, `separator`, `i`)

Splits the string `s` at matches of `separator` and returns the `i`-th of the substring between the matches.

```javascript
SplitString("Born to Be Wild", " ", 2)     // "Be"
SplitString("Born to Be Wild", "to", 1)    // " Be Wild"
SplitString("a,b,c,d", ",", 0)             // "a"
```

## <a name="startswith"></a>StartsWith(`s`, `prefix`)

Checks if the string `s` starts with a match of `prefix`.

```javascript
StartsWith("White Christmas", "White")        // true
StartsWith("White Christmas", "Christmas")    // false
```

## <a name="strfind"></a>StrFind(`s`, `sub`)

Returns the index of the first occurrence of the string `sub` in the string `s`.

If `sub` or `s` are empty then the function returns -1.

```javascript
StrFind("Hello world!", "Hello")  // 0
StrFind("Hello world!", "world")  // 6
StrFind("Hello world!", "Car")    // -1
StrFind("", "Moon")               // -1
```

## <a name="strleft"></a>StrLeft(`s`, `count`)

Returns the first `count` characters of the string `s`.

If `s` is empty then the function returns empty string.

`count` is of type integer. If `count` is negative, then the `count` number of characters will be removed from the
right-hand side of the `s` string.

```javascript
StrLeft("Daddy", 1)      // "D"
StrLeft("Daddy", 10))    // "Daddy"
StrLeft("", 3)           // ""
StrLeft("Daddy", -3)     // "Da"
StrLeft("Daddy", -8))    // ""
```

## <a name="strright"></a>StrRight(`s`, `count`)

Returns the last `count` characters of the string `s`.

If `s` is empty then the function returns empty string.

`count` is of type integer. If `count` is negative, then the `count` number of characters will be removed from the
left-hand side of the `s` string.

```javascript
StrRight("Daddy", 1)      // "y"
StrRight("Daddy", 10))    // "Daddy"
StrRight("", 3)           // ""
StrRight("Daddy", -3)     // "dy"
StrRight("Daddy", -8))    // ""
```

## <a name="strrepeat"></a>StrRepeat(`s`, `n`)

Returns a string consisting of `n` repetitions of string `s`.

```javascript
StrRepeat("Chain", 2)    // "ChainChain"
```

## <a name="substring"></a>Substring(`s`, `start`, `end`)

Returns the substring of the string `s` from `start` index, inclusive, to `end`, exclusive.

```javascript
Substring("Alfa Romeo", 5, 10)    // "Romeo"
```

## <a name="tolowercase"></a>ToLowerCase(`s`)

Converts all characters in this string `s` to lower case and returns the result obtained.

```javascript
ToLowerCase("ALPHABET")    // "alphabet"
```

## <a name="touppercase"></a>ToUpperCase(`s`)

Converts all characters in this string `s` to upper case and returns the result obtained.

```javascript
ToUpperCase("alphabet")    // "ALPHABET"
```

## <a name="trim"></a>Trim(`s`)

Returns the string without any leading and trailing whitespace of the string `s`.

```javascript
Trim("  My home.   ")    // "My home."
```
