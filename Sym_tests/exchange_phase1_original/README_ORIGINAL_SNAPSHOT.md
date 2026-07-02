# Original Snapshot

This folder is a copy of the current SKVP source files needed to build and run
the original calculation path.

No source changes were made in this folder. It is intended as the reference
code for comparison with `../exchange_phase1_modified`.

Build:

```sh
make
```

Run:

```sh
./skvp_AtomDiatom
```

Main reference file:

- `skvp_AtomDiatom.f90`: original ordered-channel calculation with `CrossSection`
  and the original `PhaseShift` path.
