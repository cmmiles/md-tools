# md-tools
A collection of Rust-based analysis tools for molecular simulations in .xtc and .trr file formats.
Full documentation can be found at [https://cmmiles.github.io/md-tools/](https://cmmiles.github.io/md-tools/).

## How to run:
`md-tools` can be run from the command line as follows:
```
md-tools order --q6 -s md-gro -t traj.xtc -o md-tools.out -start 0 --end 1000
```
where `md-tools` is the executable. The first argument after `md-tools` is the analysis to be run, options:
* `nframes` for outputting the total number of frames in the trajectory,[^note1]
* `convert` for converting into .xtc/.trr or a series of .gro/.pdb files, this can be used to shorten the trajectory,
* `order` for calculating various order parameters,
* `cluster` for calculating the largest ice-like clusters in each frame (work in progress).

### Additional arguments
* `-s` or `--sfile` \[\<.gro/.pdb\>\] structure file,
* `-t` or `--tfile` \[\<.xtc/.trr\>\] trajectory file,
* `-o` or `--output` output file,[^note2]
* `--start` \<int\> (0) analysis start time \[ps\],[^note3]
* `--end` \<int\> (+∞) analysis end time \[ps\],[^note3]
* `--stride` \<int\> (1) analysis time step \[ps\],[^note3] should be a multiple of the trajectory time step,
* `--no-pbc` disable periodic boundary conditions for analysis.
### Arguments for `md-tools order`
* `--q3` enable 3rd order Steinhardt analysis (for use with `md-tools order`),
* `--q4` enable 4th order Steinhardt analysis (for use with `md-tools order`),
* `--q6` enable 6th order Steinhardt analysis (for use with `md-tools order`).
* `--th` enable water dipole analysis (for use with `md-tools order`).
### Arguments for `md-tools cluster`
* `--n_cls` \<int\> (1) number of clusters to output for `md-tools cluster`.

[^note1]: `md-tools nframes` only accepts the arguments `-s`/`--sfile` and `-t`/`--tfile`.

[^note2]: Depending on the analysis, this may be extended, e.g. `"md-tools.out" -> "md-tools.out.q6"`.

[^note3]: Same unit as in trajectory file -- likely picoseconds.
