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
- [PadLeft](#padleft)
- [PadRight](#padright)
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

## Add(`s1`, `s2`) {#add}

Append the string `s2` at the end of the string `s1` and returns the result obtained.

```javascript
Add("Police", "man")    // "Policeman"
```

## Contains(`s`, `sub`) {#contains}

Checks if the string `s` contains a match of `sub`.

```javascript
Contains("Stormy Weather", "Sunny")      // false
Contains("Stormy Weather", "Weather")    // true
```

## EndsWith(`s`, `postfix`) {#endswith}

Checks if the string `s` ends with a match of `postfix`.

```javascript
EndsWith("White Christmas", "White")        // false
EndsWith("White Christmas", "Christmas")    // true
```

## Eq(`s1`, `s2`) {#eq}

Returns true if the string `s1` equal to `s2` and false otherwise.

```javascript
Eq("Sky", "Sky")    // true
```

## IsEmpty(`s`) {#isempty}

Returns true if the string `s` is empty or null, and false otherwise.

```javascript
IsEmpty("")         // true
IsEmpty(null)       // true
IsEmpty("Green")    // false
```

## IsNotEmpty(`s`) {#isnotempty}

Returns true if the string `s` is not empty and not null, and false otherwise.

```javascript
IsNotEmpty("")         // false
IsNotEmpty(null)       // false
IsNotEmpty("Green")    // true
```

## Length(`s`) {#length}

Returns the length of the string `s`.

```javascript
Length("Text")    // 4
```

## NotEq(`x`, `y`) {#noteq}

Returns false if the string `s1` equal to `s2` and true otherwise.

```javascript
NotEq("Sky", "Sky")    // false
```

## PadLeft(`s`, `width`, `padding`) {#padleft}

Pads the string `s` with the `padding` symbol on the left if `s` is shorter than `width`.

```javascript
PadLeft("x", 3, "-")       // "--x"
PadLeft("123", 3, "-")     // "123"
```

## PadRight(`s`, `width`, `padding`) {#padright}

Pads the string `s` with the `padding` symbol on the right if `s` is shorter than `width`.

```javascript
PadRight("x", 3, "-")       // "x--"
PadRight("123", 3, "-")     // "123"
```

## ParseFloat(`s`) {#parsefloat}

Parse `s` as a, possibly signed, real literal and return its value.

```javascript
ParseFloat("2025")      // -2025
ParseFloat("-012.78")   // -12.78
```

## ParseInt(`s`) {#parseint}

Parse `s` as a, possibly signed, integer literal and return its value.

```javascript
ParseInt("2025")    // 2025
ParseInt("-012")    // -12
```

## RegExpContains(`s`, `pattern`) {#regexpcontains}

Checks if the string `s` contains a string matches a regular expression `pattern`.

```javascript
RegExpContains("Stormy Weather", "Sunny")      // false
RegExpContains("name@gmail.com", "(\W|^)[\w.\-]{0,25}@(hotmail|gmail)\.com(\W|$)")    // true
```

## RegExpExtract(`s`, `pattern`, `n`) {#regexpextract}

Returns the `n`-th part of a string `s` that matches a regular expression `pattern`.

```javascript
RegExpExtract("Hello World!", "l+", 0)    // "ll"
RegExpExtract("Hello World!", "l+", 1)    // "l"
```

## RegExpReplace(`s`, `pattern`, `sub`) {#rexpreplace}

Returns the string after replacing a substring `sub` in string `s` according to a regular expression `pattern`.

```javascript
RegExpReplace("Hello World!", "l+", "LL")    // "HeLLo WorLLd!"
```

## ReplaceAll(`s`, `from`, `replace`) {#replaceall}

Replaces all substrings in string `s` that match `from` with `replace` and returns the result obtained.

```javascript
ReplaceAll("New York", "York", "Orleans")    // "New Orleans"
ReplaceAll("Every", "", ".")                 // ".E.v.e.r.y."
ReplaceAll("moto", "o", "")                  // "mt"
```

## SplitString(`s`, `separator`, `i`) {#splitstring}

Splits the string `s` at matches of `separator` and returns the `i`-th of the substring between the matches.

```javascript
SplitString("Born to Be Wild", " ", 2)     // "Be"
SplitString("Born to Be Wild", "to", 1)    // " Be Wild"
SplitString("a,b,c,d", ",", 0)             // "a"
```

## StartsWith(`s`, `prefix`) {#startswith}

Checks if the string `s` starts with a match of `prefix`.

```javascript
StartsWith("White Christmas", "White")        // true
StartsWith("White Christmas", "Christmas")    // false
```

## StrFind(`s`, `sub`) {#strfind}

Returns the index of the first occurrence of the string `sub` in the string `s`.

If `sub` or `s` are empty then the function returns -1.

```javascript
StrFind("Hello world!", "Hello")  // 0
StrFind("Hello world!", "world")  // 6
StrFind("Hello world!", "Car")    // -1
StrFind("", "Moon")               // -1
```

## StrLeft(`s`, `count`) {#strleft}

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

## StrRight(`s`, `count`) {#strright}

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

## StrRepeat(`s`, `n`) {#strrepeat}

Returns a string consisting of `n` repetitions of string `s`.

```javascript
StrRepeat("Chain", 2)    // "ChainChain"
```

## Substring(`s`, `start`, `end`) {#substring}

Returns the substring of the string `s` from `start` index, inclusive, to `end`, exclusive.

```javascript
Substring("Alfa Romeo", 5, 10)    // "Romeo"
```

## ToLowerCase(`s`) {#tolowercase}

Converts all characters in this string `s` to lower case and returns the result obtained.

```javascript
ToLowerCase("ALPHABET")    // "alphabet"
```

## ToUpperCase(`s`) {#touppercase}

Converts all characters in this string `s` to upper case and returns the result obtained.

```javascript
ToUpperCase("alphabet")    // "ALPHABET"
```

## Trim(`s`) {#trim}

Returns the string without any leading and trailing whitespace of the string `s`.

```javascript
Trim("  My home.   ")    // "My home."
```
