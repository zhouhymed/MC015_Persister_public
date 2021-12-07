# Cancer persistent cells

## Organization and key files

* `analysis` - the Rmarkdown files of the analyses:
    * `PC9_time_course.Rmd` - trajectory analyses of the PC9 cell line time course scRNA-seq.
    * `Survival.Rmd` - survival analysis.
* `code` - R scripts.
* `data` - R data files.
    *  `PC9_set6.merge.std.srt.Rda` - Seurat object of the down-sampled, merged and preprocessed scRNA-seq samples of different time points.
    * `PC9_set6.merge.std.session_info.v6.Rda` - session information of trajectory analysis.
    * `Survival.session_info.v6.Rda` - session information of survival analysis.
* `docs` - HTML files. You can [browse them here](https://markgene.github.io/MC015_Persister_public/index.html).

The repo is built with [workflowr][].

[workflowr]: https://github.com/jdblischak/workflowr
