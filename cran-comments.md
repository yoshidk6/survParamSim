## Test environments
* local Windows 10 install, R 3.6.2
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

Thank you for the prompt response for the initial CRAN submission.
I believe this submission addresses the comments on the initial submission.


> Thanks, if there are references describing the (theoretical background
> of) methods in your package, please add these in the Description field
> of your DESCRIPTION file in the form

There is currently no reference describing the theoretical background.


> Please replace \dontrun{} by \donttest{} or unwap the examples if they
can be executed in less than 5 sec per Rd-file.

The examples were unwrapped.


> Please replace cat() by message() or warning() in your functions (except
> for print() and summary() functions). Messages and warnings can be
> suppressed if needed.

I confirmed that all cat() were used within print() functions.
