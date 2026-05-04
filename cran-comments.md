## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on:
- Windows 11, R 4.4.x (local)

## Notes

- This is a new package submission (first release, v0.1.0).
- The `data/anomaly_benchmark.rda` file triggers an R >= 3.5.0 dependency
  note on some platforms due to the serialization format. This is expected
  and acceptable per CRAN policy.
- Vignette building requires Pandoc. The vignette uses `eval = FALSE`
  throughout so no code is executed at build time; all code is illustrative.

## Downstream dependencies

None. This is a new package with no reverse dependencies.
