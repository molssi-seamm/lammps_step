#!/usr/bin/env python
# MolSSI lammps_step:mace_mdi 1.0
"""
Optimized MACE MDI Engine for LAMMPS.

Bypasses ASE Calculator overhead by building MACE model inputs directly.
Supports GPU-accelerated neighbor lists via vesin-torch when available,
falling back to matscipy on CPU.

Usage:
    mpirun -np 1 python mace_mdi.py -mdi "-role ENGINE -name MACE -method MPI" \
        : -np 1 lmp -mdi "-role DRIVER -name LAMMPS -method MPI" -in input.dat

Authors: Paul Saxe, with assistance from Claude (Anthropic)
License: MIT
"""

import os
import sys
import time
import logging

import numpy as np
import torch
import mdi
from mpi4py import MPI
import pint

os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
logging.basicConfig(level=logging.INFO)

# Unit conversion factors that are needed
_ureg = pint.UnitRegistry()
Bohr = _ureg.Quantity(1, "bohr").to("angstrom").magnitude  # Bohr -> Å
Hartree = _ureg.Quantity(1, "hartree").to("eV").magnitude  # Hartree -> eV

# ---------------------------------------------------------------------------
# Neighbor list backends
# ---------------------------------------------------------------------------
# Try GPU-accelerated vesin-torch first, fall back to CPU matscipy
try:
    from vesin.torch import NeighborList as VesinNeighborList

    VESIN_AVAILABLE = True
    logging.info("vesin-torch available — using GPU neighbor lists")
except ImportError:
    VESIN_AVAILABLE = False
    logging.info("vesin-torch not available — using matscipy CPU neighbor lists")

if not VESIN_AVAILABLE:
    from matscipy.neighbours import neighbour_list

# Try CuEq acceleration
try:
    from mace.cli.convert_e3nn_cueq import run as run_e3nn_to_cueq

    CUEQ_AVAILABLE = True
except (ImportError, ModuleNotFoundError):
    CUEQ_AVAILABLE = False

# Try OEq acceleration
try:
    from mace.cli.convert_e3nn_oeq import run as run_e3nn_to_oeq

    OEQ_AVAILABLE = True
except (ImportError, ModuleNotFoundError):
    OEQ_AVAILABLE = False


def get_neighborhood_cpu(positions, cell, cutoff, pbc=(True, True, True)):
    """Compute neighbor list on CPU using matscipy.

    Args:
        positions: np.ndarray, shape [N, 3], in Angstroms
        cell: np.ndarray, shape [3, 3], in Angstroms
        cutoff: float, in Angstroms
        pbc: tuple of 3 bools, periodic boundary conditions

    Returns:
        edge_index: np.ndarray [2, E]
        shifts: np.ndarray [E, 3]
        unit_shifts: np.ndarray [E, 3]
    """
    cell_copy = np.array(cell, dtype=float)

    # For non-periodic directions, extend the cell so matscipy finds all neighbors
    if not all(pbc):
        identity = np.identity(3, dtype=float)
        max_pos = np.max(np.absolute(positions)) + 1
        for dim in range(3):
            if not pbc[dim]:
                cell_copy[dim, :] = max_pos * 5 * cutoff * identity[dim, :]

    sender, receiver, unit_shifts = neighbour_list(
        quantities="ijS",
        pbc=pbc,
        cell=cell_copy,
        positions=positions,
        cutoff=cutoff,
    )

    # Remove self-edges that don't cross periodic boundaries
    true_self_edge = sender == receiver
    true_self_edge &= np.all(unit_shifts == 0, axis=1)
    keep = ~true_self_edge
    sender = sender[keep]
    receiver = receiver[keep]
    unit_shifts = unit_shifts[keep]

    edge_index = np.stack((sender, receiver))
    shifts = np.dot(unit_shifts, cell_copy)

    return edge_index, shifts, unit_shifts


class MACEEngine:
    """MDI engine that drives MACE directly, bypassing ASE Calculator."""

    def __init__(
        self,
        model_path,
        device="cuda",
        default_dtype="float32",
        enable_cueq=False,
        enable_oeq=False,
    ):
        self.device = torch.device(device)
        self.dtype = torch.float32 if default_dtype == "float32" else torch.float64

        # ---- Load model ----
        model = torch.load(f=model_path, map_location=self.device, weights_only=False)

        # Convert dtype if needed
        model_dtype = next(model.parameters()).dtype
        if model_dtype != self.dtype:
            logging.warning(
                f"Model dtype {model_dtype} != requested {self.dtype}, converting."
            )
            if self.dtype == torch.float64:
                model = model.double()
            else:
                model = model.float()

        # Apply acceleration
        if enable_cueq:
            if not CUEQ_AVAILABLE:
                raise ImportError("cuequivariance not installed")
            logging.info("Converting model to CuEq for acceleration")
            model = run_e3nn_to_cueq(model, device=str(self.device)).to(self.device)
        elif enable_oeq:
            if not OEQ_AVAILABLE:
                raise ImportError("openequivariance not installed")
            logging.info("Converting model to OEq for acceleration")
            model = run_e3nn_to_oeq(model, device=str(self.device)).to(self.device)

        model.eval()
        for p in model.parameters():
            p.requires_grad_(False)

        self.model = model
        self.r_max = float(model.r_max.cpu())

        # Build z_table: maps atomic number -> one-hot index
        self.atomic_numbers = [int(z) for z in model.atomic_numbers]
        self.z_to_index = {z: i for i, z in enumerate(self.atomic_numbers)}
        self.num_species = len(self.atomic_numbers)

        # Determine head
        try:
            self.heads = list(model.heads)
        except AttributeError:
            self.heads = ["Default"]
        self.head_index = 0  # default to first head

        logging.info(
            f"Model loaded: r_max={self.r_max}, species={self.atomic_numbers}, "
            f"heads={self.heads}, dtype={self.dtype}, device={self.device}"
        )

        # ---- Pre-allocate vesin neighbor list object if available ----
        if VESIN_AVAILABLE:
            self.vesin_nl = VesinNeighborList(cutoff=self.r_max, full_list=True)
        else:
            self.vesin_nl = None

        # ---- State (set when data arrives from LAMMPS) ----
        self.natoms = None
        self.elements_np = None  # numpy, received from LAMMPS
        self.positions_np = None  # numpy, in Bohr (MDI units)
        self.cell_np = None  # numpy, in Bohr (MDI units)
        self.periodic = False  # set True when >CELL is received

        # ---- Cached tensors (allocated once, reused) ----
        self._node_attrs = None  # one-hot, on GPU
        self._batch = None  # batch index, on GPU
        self._ptr = None  # batch pointer, on GPU
        self._head = None  # head index, on GPU
        self._num_graphs = None
        self._pbc = None

        # ---- Results ----
        self.energy = None
        self.forces = None
        self.stress = None
        self._needs_calculation = True

        # ---- Timing ----
        self._n_calc = 0
        self._t_nlist = 0.0
        self._t_transfer = 0.0
        self._t_model = 0.0
        self._t_total = 0.0

    def _init_persistent_tensors(self, natoms, elements):
        """Build tensors that only change when atoms/elements change."""
        # One-hot node attributes
        indices = [self.z_to_index[int(z)] for z in elements]
        one_hot = torch.zeros(
            natoms, self.num_species, dtype=self.dtype, device=self.device
        )
        for i, idx in enumerate(indices):
            one_hot[i, idx] = 1.0
        self._node_attrs = one_hot

        # Batch tensors (single graph, always the same)
        self._batch = torch.zeros(natoms, dtype=torch.long, device=self.device)
        self._ptr = torch.tensor([0, natoms], dtype=torch.long, device=self.device)
        self._head = torch.tensor(
            [self.head_index], dtype=torch.long, device=self.device
        )
        self._num_graphs = torch.tensor(1, dtype=torch.long, device=self.device)
        # _pbc is set when periodicity is determined (on first >CELL or calculate)

    def _build_graph_vesin(self, positions_t, cell_t):
        """Build graph edges on GPU using vesin-torch."""
        # vesin expects [N, 3] positions and [3, 3] cell, both as tensors
        # It returns (edge_index [2, E], edge_vectors [E, 3])
        i, j, S, D = self.vesin_nl.compute(
            points=positions_t,
            box=cell_t,
            periodic=self.periodic,
            quantities="ijSd",
        )

        # Remove self-edges (same atom, no shift)
        self_edge = (i == j) & (S == 0).all(dim=1)
        keep = ~self_edge
        i = i[keep]
        j = j[keep]
        S = S[keep]

        edge_index = torch.stack([i, j], dim=0)
        # shifts = unit_shifts @ cell
        shifts = S.to(dtype=self.dtype) @ cell_t

        return edge_index, shifts, S.to(dtype=self.dtype)

    def _build_graph_cpu(self, positions_np, cell_np, pbc=(True, True, True)):
        """Build graph edges on CPU using matscipy, transfer to GPU."""
        edge_index_np, shifts_np, unit_shifts_np = get_neighborhood_cpu(
            positions_np, cell_np, self.r_max, pbc=pbc
        )
        edge_index = torch.tensor(edge_index_np, dtype=torch.long, device=self.device)
        shifts = torch.tensor(shifts_np, dtype=self.dtype, device=self.device)
        unit_shifts = torch.tensor(unit_shifts_np, dtype=self.dtype, device=self.device)
        return edge_index, shifts, unit_shifts

    def calculate(self):
        """Run MACE model: build graph, forward pass, extract results."""
        t_start = time.perf_counter()

        # Convert MDI Bohr -> Angstrom
        positions_ang = self.positions_np * Bohr

        if self.periodic:
            cell_ang = self.cell_np * Bohr
            pbc = (True, True, True)
            compute_stress = True
        else:
            # Non-periodic: create a large fake cell for neighbor finding
            max_pos = np.max(np.absolute(positions_ang)) + 1
            fake_size = max_pos * 5 * self.r_max
            cell_ang = np.diag([fake_size, fake_size, fake_size])
            pbc = (False, False, False)
            compute_stress = False

        # Set PBC tensor on first call
        if self._pbc is None:
            self._pbc = torch.tensor([list(pbc)], dtype=torch.bool, device=self.device)

        # ---- Neighbor list ----
        t0 = time.perf_counter()

        if self.vesin_nl is not None:
            # GPU path: transfer positions and cell, then build graph on GPU
            positions_t = torch.tensor(
                positions_ang, dtype=self.dtype, device=self.device
            )
            cell_t = torch.tensor(cell_ang, dtype=self.dtype, device=self.device)
            edge_index, shifts, unit_shifts = self._build_graph_vesin(
                positions_t, cell_t
            )
        else:
            # CPU path: build graph on CPU, then transfer
            edge_index, shifts, unit_shifts = self._build_graph_cpu(
                positions_ang, cell_ang, pbc=pbc
            )
            positions_t = torch.tensor(
                positions_ang, dtype=self.dtype, device=self.device
            )
            cell_t = torch.tensor(cell_ang, dtype=self.dtype, device=self.device)

        t1 = time.perf_counter()

        # ---- Build input dict ----
        # positions needs grad for force computation (autograd backward)
        positions_t.requires_grad_(True)

        input_dict = {
            "positions": positions_t,
            "node_attrs": self._node_attrs,
            "edge_index": edge_index,
            "shifts": shifts,
            "unit_shifts": unit_shifts,
            "cell": cell_t.unsqueeze(0),  # [1, 3, 3]
            "batch": self._batch,
            "ptr": self._ptr,
            "head": self._head,
            "num_graphs": self._num_graphs,
            "pbc": self._pbc,
        }

        t2 = time.perf_counter()

        # ---- Model forward pass ----
        out = self.model(
            input_dict,
            compute_stress=compute_stress,
            training=False,
        )

        t3 = time.perf_counter()

        # ---- Extract results, convert to MDI atomic units ----
        self.energy = out["energy"].detach().cpu().item() / Hartree
        self.forces = out["forces"].detach().cpu().to(torch.float64).numpy() / (
            Hartree / Bohr
        )

        # Stress: MACE returns [1, 3, 3] in eV/Å³, MDI expects Hartree/Bohr³
        if out.get("stress") is not None:
            self.stress = -out["stress"].detach().cpu().to(
                torch.float64
            ).numpy().reshape(3, 3) / (Hartree / Bohr**3)
        else:
            self.stress = None

        t_end = time.perf_counter()

        # ---- Timing stats ----
        self._n_calc += 1
        self._t_nlist += t1 - t0
        self._t_transfer += t2 - t1
        self._t_model += t3 - t2
        self._t_total += t_end - t_start

        if self._n_calc % 100 == 0:
            n = self._n_calc
            logging.info(
                f"Step {n}: "
                f"nlist={self._t_nlist / n * 1000:.1f}ms  "
                f"transfer={self._t_transfer / n * 1000:.1f}ms  "
                f"model={self._t_model / n * 1000:.1f}ms  "
                f"total={self._t_total / n * 1000:.1f}ms  "
                f"rate={self.natoms * n / self._t_total / 1000:.1f} katom-step/s"
            )

    # ---- MDI communication loop ----

    def run(self):
        """Main MDI engine loop."""
        mdi.MDI_Init(sys.argv[2], MPI.COMM_WORLD)

        # Register supported commands
        mdi.MDI_Register_Node("@DEFAULT")
        for cmd in [
            ">NATOMS",
            ">COORDS",
            ">CELL",
            ">ELEMENTS",
            "<ENERGY",
            "<FORCES",
            "<STRESS",
            "SCF",
            "EXIT",
        ]:
            mdi.MDI_Register_Command("@DEFAULT", cmd)

        comm = mdi.MDI_Accept_Communicator()
        logging.info("MDI connection established")

        while True:
            command = mdi.MDI_Recv_Command(comm)
            logging.debug(f"MDI command: {command}")

            if command == "EXIT":
                break

            elif command == ">NATOMS":
                self.natoms = mdi.MDI_Recv(1, mdi.MDI_INT, comm)

            elif command == ">ELEMENTS":
                elements = mdi.MDI_Recv(self.natoms, mdi.MDI_INT, comm)
                self.elements_np = np.array(elements, dtype=np.int64)
                self._init_persistent_tensors(self.natoms, self.elements_np)
                logging.info(
                    f"Received {self.natoms} atoms, "
                    f"elements: {sorted(set(self.elements_np.tolist()))}"
                )

            elif command == ">CELL":
                cell = mdi.MDI_Recv(9, mdi.MDI_DOUBLE, comm)
                self.cell_np = np.array(cell, dtype=np.float64).reshape(3, 3)
                if not self.periodic:
                    self.periodic = True
                    self._pbc = torch.tensor(
                        [[True, True, True]], dtype=torch.bool, device=self.device
                    )
                    logging.info("Periodic system detected")
                self._needs_calculation = True

            elif command == ">COORDS":
                coords = mdi.MDI_Recv(3 * self.natoms, mdi.MDI_DOUBLE, comm)
                self.positions_np = np.array(coords, dtype=np.float64).reshape(
                    self.natoms, 3
                )
                self._needs_calculation = True
                if self._n_calc < 2:
                    logging.debug(f"First 3 positions (Bohr): {self.positions_np[:3]}")
                    self._needs_calculation = True

            elif command == "<ENERGY":
                if self._needs_calculation:
                    self.calculate()
                    self._needs_calculation = False
                mdi.MDI_Send(self.energy, 1, mdi.MDI_DOUBLE, comm)

            elif command == "<FORCES":
                if self._needs_calculation:
                    self.calculate()
                    self._needs_calculation = False
                mdi.MDI_Send(
                    self.forces.flatten(), 3 * self.natoms, mdi.MDI_DOUBLE, comm
                )

            elif command == "<STRESS":
                if self._needs_calculation:
                    self.calculate()
                    self._needs_calculation = False
                if self.stress is not None:
                    if self._n_calc < 2:
                        logging.debug(f"MDI sending {self.stress=}")
                    mdi.MDI_Send(self.stress.flatten(), 9, mdi.MDI_DOUBLE, comm)
                else:
                    # Send zeros if stress wasn't computed
                    if self._n_calc < 2:
                        logging.debug("MDI sending zeroes for stress")
                    mdi.MDI_Send(np.zeros(9), 9, mdi.MDI_DOUBLE, comm)

            elif command == "SCF":
                self.calculate()
                self._needs_calculation = False

            else:
                print(f"Error: unhandled MDI command {command}!", file=sys.stderr)
                sys.exit(1)

        logging.info(
            f"Engine finished. {self._n_calc} calculations, "
            f"avg {self._t_total / max(self._n_calc, 1) * 1000:.1f} ms/step"
        )

        # Clean up PyTorch before MPI tears down the process
        import gc

        torch.cuda.synchronize()
        del self.model
        del self._node_attrs
        gc.collect()
        torch.cuda.empty_cache()


if __name__ == "__main__":
    model_path = os.environ.get("SEAMM_FF")
    if model_path is None:
        print("Error: SEAMM_FF environment variable not set", file=sys.stderr)
        sys.exit(1)

    engine = MACEEngine(
        model_path,
        device="cuda:0",
        default_dtype="float32",
        enable_cueq=True,
    )
    engine.run()
