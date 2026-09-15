Machine-learned forcefields: one engine, chosen by configuration
================================================================

Notes for the campaign that made a machine-learned forcefield run reach its
engine without the plug-in knowing which engine that is, and that settled
``xnns`` as the single mechanism for them.

The starting point was a job that hung for four hours and died on a
``torch.load``. Getting from there to a working MLFF run turned up a chain of
faults, most of which were only reachable because the one before it had been
fixed. The notes below record the design that came out of it, and the reasoning
that is not visible in the code.

Contents:

.. toctree::
   :glob:
   :maxdepth: 2

   NOTES_*
