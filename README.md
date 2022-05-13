# SCONCE2

This directory is for the program SCONCE2.

## Dependencies
SCONCE2 is written in C++11, and requires
- GNU make (tested on v4.1)
- g++ (tested on 7.5.0)
- BOOST (tested on v1.66.0)
- GSL (tested on v2.4)

Additional [R scripts](scripts/) require
- R (tested on v4.2.0)
- ape (tested on v5.6-2)
- cowplot (tested on v1.1.1)
- ggplot2 (tested on v3.3.5)
- ggtree (tested on v3.4.0)
- grid (tested on v0.5-1)
- gtools (tested on v3.9.2)
- phangorn (tested on v2.8.1)
- plyr (tested on v1.8.7)
- reshape2 (tested on v1.4.4)
- scales (tested on v1.2.0)
- stringr (tested on v1.4.0)

SCONCE2 was developed and tested on Ubuntu 18.04.

## Installation instructions
1. Clone this repo:
```
git clone git@github.com:NielsenBerkeleyLab/sconce2.git
```
2. Run `make`. This will build intermediates into the `build/` directory and create an executable named `sconce2`.

