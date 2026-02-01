# How to contribute to OpenFMS

Thanks for contributing to OpenFMS!💜
Here's a couple of guidelines that you should keep in mind.

## Setup your environment

We use tools such as autoformatter ([fprettify](https://github.com/fortran-lang/fprettify)) and linter ([fortitude](https://github.com/PlasmaFAIR/fortitude/))
to keep our code nice and tidy. Instead of having the developer to
install these and run them manually, we use a tool called `prek`
which does this automatically on every commit.

To install `prek` and the `configargparse` dependency (needed for `fprettify`) run:

```console
pip install --user prek configargparse
```

Then you must run the following in the repository to install prek's [pre-commit hook](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks):

```
prek install
```

From now on, prek will execute hooks defined in `.pre-commit-config.yaml` before every commit.
You can also run the hooks manually at any point:

```
❯ prek run -a
don't commit to branch.................................Passed
check yaml.............................................Passed
check for added large files............................Passed
check for merge conflicts..............................Passed
check that executables have shebangs...................Passed
fortitude..............................................Passed
Check module usage.....................................Passed
Autoformat code........................................Passed
```

To run only autoformatter

```
❯ prek run -a format
Autoformat code........................................Passed
```

See [prek documentation](https://prek.j178.dev/) for more details.


## Code style

Here's a quick summary of our code style that we try to adhere to (especially in new code):

 - use `implicit none` everywhere, no exceptions! This is enforced by the compiler in CI.
 - each function/subroutine argument should have the intent attribute.
 - use `real(DefReal)` for real numbers, the `DefReal` constant indicating the kind (precision) is defined in `GlobalModule.f90`
 - variables and subroutines in modules should be private by default, use `private` attribute,e.g.

```fortran
module mod_my_module
   private
   integer :: public_var
   integer :: private_var
   public :: public_var
   public :: public_subroutine
contains
   subroutine public_subroutine
      ...
   end subroutine public_subroutine
end module mod_my_module
```

### Code formatting

Code formatting is enforced using `fprettify` formatter,
and you don't need to think about it for the most part.
Here's a quick summary of our formatting style, as it is defined in `.fprettify.rc` config file.

 - use 3 spaces for indentation, no tabs.
 - use capital letters only for defined constants (e.g. those from module `mod_const`). We use lower case for everything else.
 - use C-style relational operators (`< > == /=`) instead of the old FORTRAN style (`.gt. .lt.`)
 - comments should start at the same indentation level as the code they are commenting.
    - use an exclamation mark to start a comment

## Submitting code changes

Last but not least, to get your code merged to the main repository, please open a Pull Request (PR) on Github.
If you're not familiar with Pull Requests, take a look [here](https://guides.github.com/activities/hello-world/#pr).

It's super easy, barely an inconvenience! (assuming basic familiarity with Git)

### Inspecting Git history

To ignore bulk insignificant changes in git blame history (e.g. whitespace changes),
configure `git blame` to ignore commits specified in `.git-blame-ignore-revs` file.

```sh
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

This is done automatically in GitHub UI.
