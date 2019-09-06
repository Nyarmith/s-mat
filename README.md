## s-mat: a static matrix library
Personal matrix library that tries to do everything at compile-time if possible.

To use, include `include/s-mat.hh`, or with cmake `target_link_library(<an_exe> s-mat)`, define `ENABLE_TESTING` in cmake to build tests

TODO:
- make matrices iterable in row-major
- create notation for applying functions element-wise (or just a separate fn)
- when complete enough, add actual usage instructions in this readme
