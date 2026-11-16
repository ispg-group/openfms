# How to contribute to OpenFMS

Thanks for contributing to OpenFMS!💜
Here's a couple of guidelines to keep in mind.

## Setup your environment

During development we use tools such as autoformatter ([fprettify](https://github.com/fortran-lang/fprettify)) and linter ([fortitude](https://github.com/PlasmaFAIR/fortitude/))
to keep our code nice and tidy. Instead of having the developer to
install these and run them manually, we use a program called `prek`
which does this automatically on every commit.

1. [Install prek](https://prek.j178.dev/installation/)

2. Run `prek install` in your local OpenFMS repository to install prek's [pre-commit hook](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks):

From now on, prek will execute hooks defined in `prek.toml` before every commit.
You can also run the hooks manually at any point:

```
❯ prek run --all-files
don't commit to branch.................................Passed
check yaml.............................................Passed
check for added large files............................Passed
check for merge conflicts..............................Passed
check that executables have shebangs...................Passed
fortitude..............................................Passed
Check module usage.....................................Passed
Autoformat code........................................Passed
```

To run only autoformatter (`--all-files` abbreviated as `-a`)

```
❯ prek run -a format
Autoformat code........................................Passed
```

If any of the checks fail, the commit is aborted. Often times, the violations
are fixed automatically (e.g. formatter will autoformat the code), so simply
re-running `git commit` together with the changed files is enough.
Sometimes, manual intervention is necessary, for example if `fortitude` catches
that you forgot to use `implicit none`.

```console
❯ prek run -a fortitude
fortitude................................................................Failed
- hook id: fortitude
- exit code: 1

  src/abin.F90:19:1: C001 program missing 'implicit none'
     |
  17 | !  You should have received a copy of the GNU General Public License
  18 | !  along with this program in the file LICENSE. If not, see <http://www.gnu.org/licenses/>.
  19 | program abin
     | ^^^^^^^^^^^^ C001
  20 |    use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
  21 |    use mod_const, only: DP, AUtoFS
  22 |    use mod_arrays
     |

  fortitude: 47 files scanned.
  Number of errors: 1
```


> [!TIP]
> In rare circumstances, you might want to skip the checks and make a commit even if they are failing,
> such as when you want to quickly save your work.
> Use the `--no-verify` (or `-n`) option to skip the pre-commit hooks.
> `git commit --no-verify`
> Note that the pre-commit hooks are nevertheless always enforced on GitHub.

See [prek documentation](https://prek.j178.dev/) for more details.


## Code style

Here's a quick summary of our code style that we try to adhere to (especially in new code):

 - use `implicit none` everywhere, no exceptions! This is enforced by the compiler in CI.
 - each function/subroutine dummy argument should have the intent attribute.
 - use `real(DefReal)` for real numbers, `integer(DefInt)` for integers and complex(DefComp) for complex numbers.
   The `DefXXX` constants indicating the kind (precision) are defined in `GlobalModule`.
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
Here's a quick summary of our formatting style, defined in `.fprettify.rc` config file.

 - use 3 spaces for indentation, no tabs.
 - use C-style relational operators (`< > == /=`) instead of the old FORTRAN style (`.gt. .lt.`)
 - comments should start at the same indentation level as the code they are commenting.
    - use an exclamation mark to start a comment

## Submitting code changes

Last but not least, to get your code merged to the main repository, please open a Pull Request (PR) on Github.
If you're not familiar with Pull Requests, [take a look here](https://guides.github.com/activities/hello-world/#pr).

### Inspecting Git history

To ignore bulk insignificant changes in git blame history (e.g. whitespace changes),
configure `git blame` to ignore commits specified in `.git-blame-ignore-revs` file.

```sh
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

This is done automatically in GitHub UI.

## Building the developer documentation

Building the docs is straightforward:

1. [Install FORD with your favorite python package manager](https://forddocs.readthedocs.io/en/stable/index.html#sec-installation)
2. Execute `ford ford_config.md`

That's it.
