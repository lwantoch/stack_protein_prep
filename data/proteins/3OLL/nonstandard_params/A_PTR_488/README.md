# Pre-existing Parameters Found — PTR  |  3OLL / A_PTR_488

Parameters for **PTR** are already available in **phosaa14SB**.

phosphotyrosine (−2)

## Loading in tleap

```tcl
source leaprc.protein.ff14SB
source leaprc.phosaa14SB
# PTR is now available as a recognised residue
```

## No RESP derivation needed

Because pre-built, peer-reviewed parameters exist, you do not need to run
Gaussian for this residue.  If you want to derive independent RESP parameters
(e.g. for publication or to ensure ff14SB consistency), a scaffold is provided
in the ``resp/`` subdirectory.
