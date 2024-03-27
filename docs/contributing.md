# Contribute

Contributions, either with [issues](https://github.com/pierre-24/stdlite/issues) or [pull requests](https://github.com/pierre-24/stdlite/pulls) are welcomed.

## Install

If you want to contribute, this is the usual deal: 
start by [forking](https://guides.github.com/activities/forking/), then clone your fork and [install the dependencies using Meson](build.md).

I personally use [Clion](https://www.jetbrains.com/clion/) for the development, as it is integrated with the Meson build system, but others might do.

## Tips to contribute

+ A good place to start is the [list of issues](https://github.com/pierre-24/stdlite/issues).
  In fact, it is easier if you start by filling an issue, and if you want to work on it, says so there, so that everyone knows that the issue is handled.

+ Don't forget to work on a separate branch.
  Since this project follow the [git flow](http://nvie.com/posts/a-successful-git-branching-model/), you should base your branch on `dev`, not work in it directly:

    ```bash
    git checkout -b new_branch origin/dev
    ```

+ Pull requests should be unitary, and include unit test(s) and documentation if needed.
  The test suite must succeed for the merge request to be accepted.

+ If you want to see and edit the doc, you can run the `mkdocs` webserver:

    ```bash
    pip install mkdocs mkdocs-material
    mkdocs serve
    ```
  
    Code documentation is generated via [doxide](https://doxide.org/), with:
    
    ```bash
    doxide build
    ```