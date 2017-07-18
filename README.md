# repovar: report minority variants

This project is a spin-off of the minority variants software [MinVar].

`repovar` takes care of the creation of a report (in markdown first, then
converted with pandoc into pdf). This will make the configuration of the
report easier also for the third user.

#### Simple usage

After a run with MinVar one has (among the others) csv files in the output
directory:
- `annotated_DRM.csv` and
- `subtype_evidence.csv`.

In the output directory, one runs

    repovar -m annotated_DRM.csv -s subtype_evidence.csv

and a markdown file `report.md` is created. Here below the begin of this file
is reported.

    Drug resistance mutations detected by NGS sequencing
    ====================================================

    Subtype inference with blast
    ----------------------------
    |     subtype     | support [%] |
    |:---------------:|:-----------:|
    | CONSENSUS_02_AG |     91      |
    |     A1.anc      |      4      |
    |   CON_OF_CONS   |      1      |
    |     Mgroup      |      1      |
    |CONSENSUS_06_CPX |      1      |


    Parsing mutations
    -----------------

    The list of mutations was downloaded from HIVdb and includes:


## Make a pdf report

For the creation of a pdf report [pandoc] and pdflatex must be installed.
pandoc is available in many ways. pdflatex is available via the
TexLive installation for your system, see the instructions on [TexLive] page.

[MinVar]: https://git.io/MinVar
[pandoc]: https://www.pandoc.org
[TexLive]: https://www.tug.org/texlive/
