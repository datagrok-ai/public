<!-- TITLE: Markdown -->
<!-- SUBTITLE: -->

# Markdown

[Markdown](http://daringfireball.net/projects/markdown) is a way to style text on the web. 
You control the display of the document; formatting words as bold or italic, adding images, 
and creating lists are just a few of the things we can do with Markdown. Mostly, Markdown 
is just regular text with a few non-alphabetic characters thrown in, like `#` or `*`.

## Syntax guide

Here’s an overview of Markdown syntax that you can use.

### Headers

    # This is an <h1> tag
    ## This is an <h2> tag
    ###### This is an <h6> tag
    

### Emphasis

    *This text will be italic*
    _This will also be italic_
    
    **This text will be bold**
    __This will also be bold__
    
    _You **can** combine them_
    
### Lists

#### Unordered

    * Item 1
    * Item 2
      * Item 2a
      * Item 2b
    
#### Ordered

    1. Item 1
    1. Item 2
    1. Item 3
       1. Item 3a
       1. Item 3b
    
### Images

    ![GitHub Logo](/images/logo.png)
    Format: ![Alt Text](url)
    
### Links

    http://github.com - automatic!
    [GitHub](http://github.com)
    
### Blockquotes

    As Kanye West said:
    
    > We're living the future so
    > the present is our past.
    
### Inline code

    I think you should use an
    `<addr>` element here instead.
    
### Syntax highlighting

Here’s an example of how you can use syntax highlighting with 
[GitHub Flavored Markdown](https://help.github.com/articles/basic-writing-and-formatting-syntax):

    ```javascript
    function fancyAlert(arg) {
      if(arg) {
        $.facebox({div:'#foo'})
      }
    }
    ```
    
You can also simply indent your code by four spaces:

        function fancyAlert(arg) {
          if(arg) {
            $.facebox({div:'#foo'})
          }
        }
    
Here’s an example of Python code without syntax highlighting:

    def foo():
        if not bar:
            return True
    

### Task Lists

    - [x] @mentions, #refs, [links](), **formatting**, and <del>tags</del> supported
    - [x] list syntax required (any unordered or ordered list supported)
    - [x] this is a complete item
    - [ ] this is an incomplete item
    
If you include a task list in the first comment of an Issue, you will get a handy progress 
indicator in your issue list. It also works in Pull Requests!

### Tables

You can create tables by assembling a list of words and dividing them with 
hyphens `-` (for the first row), and then separating each column with a pipe `|`:

    | First Header                | Second Header                |
    |-----------------------------|------------------------------|
    | Content from cell 1         | Content from cell 2          |
    | Content in the first column | Content in the second column |
    
    
### Datagrok-specific Markup

A powerful mechanism that allows to embed platform-specific visual elements right into 
the HTML document. Simply embed the expression like that: `#{expression}`, see more 
information [there](markup.md).


See also:
  * [Mastering Markdown](https://guides.github.com/features/mastering-markdown)
  * [GitHub Flavored Markdown](https://help.github.com/articles/basic-writing-and-formatting-syntax)
