# Create Metagene intervals
"metagene" intervals associated with a given interval (input, in bed format) is the set of intervals positioned with respect to the input interval. This allows for metagene plots to be made.

## Features
- Strand-aware
- normalize different input interval length
- option to center on the middle of interval
- option to fold with respect to boundary

## Usage

```sh
python ~/wzlib/pyutils/wzmetagene.py interval_bed_file
```

## Examples
#### Sample 3 intervals internal of the input interval evenly.
```sh
$ echo -e "chr1\t10000\t10200" | python wzmetagene.py - -n 3 -m 0
chr1    10000   10066   0       0-33%   0       chr1    10000   10200
chr1    10066   10133   1       33-66%  0       chr1    10000   10200
chr1    10133   10200   2       66-100% 0       chr1    10000   10200
```

#### Sample 3 intervals internal and 2 intervals flanking the input interval
The flanking sequence is 1kb in length and specified by `-f`.
```sh
$ echo -e "chr1\t10000\t10200" | python wzmetagene.py - -n 3 -m 2 -f 1000
chr1    9500    10000   -1      (-500)-(0)      -1      chr1    10000   10200
chr1    9000    9500    -2      (-1000)-(-500)  -1      chr1    10000   10200
chr1    10000   10066   0       0-33%   0       chr1    10000   10200
chr1    10066   10133   1       33-66%  0       chr1    10000   10200
chr1    10133   10200   2       66-100% 0       chr1    10000   10200
chr1    10200   10700   3       0-500   1       chr1    10000   10200
chr1    10700   11200   4       500-1000        1       chr1    10000   10200
```

#### Input interval has a strand given in an extra column
The strand column is given in `-s`.
```sh
$ echo -e "chr1\t10000\t10200\t-" | python wzmetagene.py - -n 3 -m 2 -f 1000 -s 4
chr1    10200   10700   -1      (-500)-(0)      -1      chr1    10000   10200   -
chr1    10700   11200   -2      (-1000)-(-500)  -1      chr1    10000   10200   -
chr1    10000   10066   2       66-100% 0       chr1    10000   10200   -
chr1    10066   10133   1       33-66%  0       chr1    10000   10200   -
chr1    10133   10200   0       0-33%   0       chr1    10000   10200   -
chr1    9500    10000   3       0-500   1       chr1    10000   10200   -
chr1    9000    9500    4       500-1000        1       chr1    10000   10200   -
```

#### Report the middle coordinate (1 in length) as opposed to full output interval
This is done by `--middle`
```sh
$ echo -e "chr1\t10000\t10200\t-" | python wzmetagene.py - -n 3 -m 0 --middle
chr1    10032   10033   0       0-33%   0       chr1    10000   10200   -
chr1    10099   10100   1       33-66%  0       chr1    10000   10200   -
chr1    10165   10166   2       66-100% 0       chr1    10000   10200   -
```

#### Collapse the input interval and make flanking sequence from the middle
This is achieved by `--collapse`. Note it's different from `--middle` which only change the output intervals.
`-F` specifies the step instead of specifying the whole flanking length as by `-f`.
```sh
$ echo -e "chr1\t10000\t10200" | python wzmetagene.py - -n 3 -F 400 -m 3 --collapse
chr1    9700    10100   -1      (-400)-(0)      -1      chr1    10000   10200
chr1    9300    9700    -2      (-800)-(-400)   -1      chr1    10000   10200
chr1    8900    9300    -3      (-1200)-(-800)  -1      chr1    10000   10200
chr1    10100   10500   0       0-400   1       chr1    10000   10200
chr1    10500   10900   1       400-800 1       chr1    10000   10200
chr1    10900   11300   2       800-1200        1       chr1    10000   10200
```

##### Like above but also have an output interval centered on the middle of the input interval
This is achieved by `--collapsecenter`.
```sh
$ echo -e "chr1\t10000\t10200" | python wzmetagene.py - -n 3 -F 400 -m 3 --collapsecenter
chr1    9500    9900    -1      (-400)-(0)      -1      chr1    10000   10200
chr1    9100    9500    -2      (-800)-(-400)   -1      chr1    10000   10200
chr1    8700    9100    -3      (-1200)-(-800)  -1      chr1    10000   10200
chr1    9900    10300   0       0-100%  0       chr1    10000   10200
chr1    10300   10700   1       0-400   1       chr1    10000   10200
chr1    10700   11100   2       400-800 1       chr1    10000   10200
chr1    11100   11500   3       800-1200        1       chr1    10000   10200
```

##### Create a symmetric folding behavior with respect to the boundary
This is done by `--fold`
```sh
$ echo -e "chr1\t10000\t10200" | python wzmetagene.py - -n 4 -f1000 -m 4 --fold
chr1    9750    10000   -1      (-250)-(0)      -1      chr1    10000   10200
chr1    9500    9750    -2      (-500)-(-250)   -1      chr1    10000   10200
chr1    9250    9500    -3      (-750)-(-500)   -1      chr1    10000   10200
chr1    9000    9250    -4      (-1000)-(-750)  -1      chr1    10000   10200
chr1    10000   10050   0       0-25%   0       chr1    10000   10200
chr1    10050   10100   1       25-50%  0       chr1    10000   10200
chr1    10100   10150   1       25-50%  0       chr1    10000   10200
chr1    10150   10200   0       0-25%   0       chr1    10000   10200
chr1    10200   10450   -1      (-250)-(0)      -1      chr1    10000   10200
chr1    10450   10700   -2      (-500)-(-250)   -1      chr1    10000   10200
chr1    10700   10950   -3      (-750)-(-500)   -1      chr1    10000   10200
chr1    10950   11200   -4      (-1000)-(-750)  -1      chr1    10000   10200
```