Forcefields as model chemistries: retiring the Forcefield step
==============================================================

Notes for the campaign that makes the LAMMPS step take its interatomic potential
from the Model Chemistry step -- classical forcefields, ReaxFF, OpenKIM
potentials and machine-learned forcefields alike -- and retires the Forcefield
step as a step in flowcharts.

It started from a question about why production MLFF MD and labelling ran from
hand-written scripts rather than through SEAMM (2026-10-11). Part of the answer
was that a fresh MLFF checkpoint, with its per-run engine options, is awkward to
reach through the Forcefield step. Paul's conclusion: the Forcefield step does
nothing a model chemistry cannot do.

Contents:

.. toctree::
   :glob:
   :maxdepth: 2

   NOTES_*
