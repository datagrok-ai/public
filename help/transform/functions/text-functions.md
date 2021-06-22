<!-- TITLE: Text functions -->
<!-- SUBTITLE: -->

# Text functions

## RegExpExtract(`s`, `pattern`, `n`)

Returns the `n`-th part of a string `s` that matches a regular expression `pattern`.

```javascript
/* Example of using the RegExpExtract */
RegExpExtract('Hello World!', 'l+', 1)    // Returns `ll`
RegExpExtract('Hello World!', 'l+', 2)    // Returns `l`
```

## RegExpReplace(`s`, `pattern`, `sub`)

Returns the string after replacing a substring `sub` in string `s` according to a regular expression `pattern`.

```javascript
/* Example of using the RegExpReplace */
RegExpReplace('Hello World!', 'l+', 'LL')    // Returns `HeLLo WorLLd!`
```

## StrFind(`s`, `sub`)

Returns the index of the first occurrence of the string `sub` in the string `s`.

Indexing characters in a string starts at one. If `sub` or `s` are empty then the function returns 0.

```javascript
/* Examples of using the StrFind */
StrFind('Hello world!', 'Hello')  // Returns 1
StrFind('Hello world!', 'world')  // Returns 7
StrFind('Hello world!', 'Car')    // Returns 0
StrFind('', 'Moon')               // Returns 0
```

## StrLeft(`s`, `count`)

Returns the first `count` characters of the string `s`.

Indexing characters in a string starts at one. If `s` is empty then the function returns empty string.

`count` can be either real or integer type. If `count` is of type real only the integer part is used. If `count` is negative, then the `count` number of characters will be taken from the right-hand side of the `s` string. If `count` is greater than the length of `s`, the whole string is returned.

```javascript
/* Examples of using the StrLeft */
StrLeft('Daddy', 1)      // Returns 'D'
StrLeft('Daddy', 5.7)    // Returns 'Daddy'
StrLeft('Daddy', 10))    // Returns 'Daddy'
StrLeft('Daddy', 0.1)    // Returns ''
StrLeft('', 3)           // Returns ''
StrLeft('Daddy', -3)     // Returns 'ddy'
StrLeft('Daddy', -8))    // Returns 'Daddy'
```

## StrRight(`s`, `count`)

Returns the last `count` characters of the string `s`.

Indexing characters in a string starts at one. If `s` is empty then the function returns empty string.

`count` can be either real or integer type. If `count` is of type real only the integer part is used. If `count` is negative, then the `count` number of characters will be taken from the left-hand side of the `s` string. If `count` is greater than the length of `s`, the whole string is returned.

```javascript
/* Examples of using the StrRight */
StrRight('Daddy', 1)      // Returns 'y'
StrRight('Daddy', 5.7)    // Returns 'Daddy'
StrRight('Daddy', 10))    // Returns 'Daddy'
StrRight('Daddy', 0.1)    // Returns ''
StrRight('', 3)           // Returns ''
StrRight('Daddy', -3)     // Returns 'Dad'
StrRight('Daddy', -8))    // Returns 'Daddy'
```

## StrRepeat(`s`, `num`)

Returns a string consisting of `num` repetitions of string `s`.

```javascript
/* Example of using the StrRepeat */
StrRepeat('Chain', 2)    // Returns 'ChainChain'
```

## Trim(`s`)

Returns the string without any leading and trailing whitespace of the string `s`.

```javascript
/* Example of using the Trim */
Trim('  My home.   ')    // Returns 'My home.'

```
