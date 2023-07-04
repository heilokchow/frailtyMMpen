This is an update for this package

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## Previous cran-comments

## Resubmission

This is a resubmission. In this version I have:

* Explain all acronyms in the description tex and remove the single quote for these acronym

## Resubmission

This is a resubmission. In this version I have:

* Edit the description file
* Reduce the length of the title to less than 65 characters.
* Change MM -> 'MM' in description.
* Add reference in authors (year) <doi:...> in the description field of DESCRIPTION file.
* Add \value to the following Rd files which previously miss \value tag:
      cluster.Rd: \value
      coef.fpen.Rd: \value
      event.Rd: \value
      print.fmm.Rd: \value
      print.fpen.Rd: \value
* Replace \dontrun{} by \donttest{}.

## Test environments

* Local:
  - Windows 10, R 4.2.2 (x86_64-w64-mingw32/x64 (64-bit))
* win-builder:
  - Windows Server 2022 (release and devel)
    - 1 note since this is the first submission from my email address and the doi and citations are correct:
      ```
      * checking CRAN incoming feasibility ... NOTE
      Maintainer: 'Yunpeng Zhou <u3514104@connect.hku.hk>'
      
      New submission
      
      Possibly misspelled words in DESCRIPTION:
      MCP (16:245)
      Minorize (16:46)
      Xu (17:12, 18:12)
      Zhou (17:19, 18:19)
      minimax (16:220)
      
      ```
* R-hub builder (https://builder.r-hub.io)
  - Windows Server 2022 R-devel, 64 bit
    - 1 note same as win-builder and other notes are unrelated to my package:
      ```
      * checking CRAN incoming feasibility ... NOTE
      Maintainer: 'Yunpeng Zhou <u3514104@connect.hku.hk>'
      
      New submission
      * checking HTML version of manual ... NOTE
      Skipping checking math rendering: package 'V8' unavailable
      * checking for detritus in the temp directory ... NOTE
      Found the following files/directories:
        'lastMiKTeXException'
      ```
* Github Actions
  - Ubuntu-20.04

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.
