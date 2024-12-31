# libcirclewire
Library to process 2D wire loops consisting of straight lines and circular arcs.  Can also process extrusions and revolutions of 2D wire loops.

# Dependancies

1. [libpolyhedra](https://github.com/maurerpe/libpolyhedra)
2. GNU Autotools and a functional C compilier.

# System requirements

POSIX compliant OS (such as Linux or FreeBSD) or Windows.  Really, any system that can run the dependancies.

# Installation

```
autoconf -i
./configure
make
sudo make install
```

Run `./configure --help` for configuration options.
