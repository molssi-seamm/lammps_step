#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The executor's result is logged without the contents of the returned files."""

import pprint

from lammps_step.lammps import _summarize_run_result


def big_result(n=5_000_000):
    return {
        "stdout": "LAMMPS (2 Aug 2023)\n" * (n // 20),
        "stderr": "",
        "files": ["com_velocities.trj", "log.cite", "broken.trj"],
        "com_velocities.trj": {
            "exception": None,
            "data": "1 0.1 0.2 0.3\n" * (n // 14),
        },
        "log.cite": {"exception": None, "data": b"citation"},
        "broken.trj": {"exception": "IOError", "data": None},
        "listing": "a\\nb",
    }


def test_files_are_summarized_by_size():
    result = big_result()
    summary = _summarize_run_result(result)
    assert summary["com_velocities.trj"] == {
        "size": len(result["com_velocities.trj"]["data"])
    }
    assert summary["log.cite"] == {"size": 8}
    assert summary["broken.trj"] == {"size": None, "exception": "IOError"}
    assert summary["stdout"].endswith(f"... ({len(result['stdout'])} characters)")
    assert summary["files"] == result["files"] and summary["stderr"] == ""
    # The result itself is untouched
    assert len(result["com_velocities.trj"]["data"]) > 1_000_000
    assert len(pprint.pformat(summary)) < 5000


def test_files_given_by_path_are_summarized(tmp_path):
    path = tmp_path / "com_velocities.trj"
    path.write_text("1 0.1 0.2 0.3\n" * 1000)
    result = {
        "files": ["com_velocities.trj"],
        "com_velocities.trj": {"path": path, "exception": None},
    }
    summary = _summarize_run_result(result)
    assert summary["com_velocities.trj"] == {
        "path": str(path),
        "size": path.stat().st_size,
    }
