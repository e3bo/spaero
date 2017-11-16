
Contributions may take the form of GitHub issues or pull
requests. Issues can be used to start discussions, request
clarification in the documentation, or suggest features. However,
another key type of issue is a bug report. To submit a great bug
report, please include a small reproducible example of the
bug. [This Stack Overflow page](https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example)
describes many techniques for doing this.

Pull requests can be used to submit changes to the project's files for
consideration of inclusion in the project. These changes can be as
small as fixing a typo or clarifying a sentence in the documentation
or as large as the implementation of a major new feature. This project
uses the typical GitHub pull request workflow. The Pro Git ebook
contains a nice introduction to the steps involved
[here](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project).

Pull requests will be evaluated according to the following criteria:

1. __relevance__. The goal of the changes must be clearly stated and
   consistent with spaero's general goals.
2. __coherence__. The changes in the pull request must all be related
   to the goal of the pull request.
4. __trustworthiness__. Code should be written in a modular pattern
   that is easy to read and check with tests. Within reason, all lines of
   code should be exercised by tests, as evaluated by the covr
   package. New tests should be added to the existing testing
   structure, which is based on the guidelines in the
   [R Packages book](http://r-pkgs.had.co.nz/tests.html).
4. __documentation__. Any new parameters, functions, or function
   outputs made available to the user must be documented. The
   documentation for the project uses Roxygen for the manual and
   rmarkdown for vignettes as described in the
   [R Packages book](http://r-pkgs.had.co.nz/man.html).
3. __style__. Code style should be consistent with the default
   style that is enforced by the lintr package.

It also often helps the maintainers evaluate pull requests when
individual commits in the pull request are atomic and have messages
that explain why changes were made.
