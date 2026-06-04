"""Tests for gaussian_hpc.SlurmConfig — CESGA FinisTerrae III CPU-only defaults."""
from __future__ import annotations

from stack_protein_preparation.gaussian_hpc import SlurmConfig


class TestSlurmConfigDefaults:
    def test_partition_is_short(self) -> None:
        assert SlurmConfig().partition == "short"

    def test_time_is_six_hours(self) -> None:
        assert SlurmConfig().time == "06:00:00"

    def test_cpus_is_16(self) -> None:
        assert SlurmConfig().cpus == 16

    def test_gpus_is_zero(self) -> None:
        assert SlurmConfig().gpus == 0

    def test_mem_default(self) -> None:
        assert SlurmConfig().mem == "32G"


class TestSlurmConfigResourceArgs:
    def test_partition_flag_present(self) -> None:
        assert "--partition=short" in SlurmConfig().resource_args()

    def test_cpus_per_task_flag(self) -> None:
        assert "--cpus-per-task=16" in SlurmConfig().resource_args()

    def test_mem_is_total_not_per_cpu(self) -> None:
        args = SlurmConfig().resource_args()
        assert any(a.startswith("--mem=") for a in args)
        assert not any(a.startswith("--mem-per-cpu=") for a in args)

    def test_no_gres_by_default(self) -> None:
        args = SlurmConfig().resource_args()
        assert not any("gres" in a for a in args)

    def test_gres_a100_when_gpus_set(self) -> None:
        args = SlurmConfig(gpus=1).resource_args()
        assert "--gres=gpu:a100:1" in args

    def test_no_requeue_flag(self) -> None:
        assert "--requeue" not in SlurmConfig().resource_args()

    def test_nodes_and_ntasks_are_one(self) -> None:
        args = SlurmConfig().resource_args()
        assert "--nodes=1" in args
        assert "--ntasks=1" in args

    def test_custom_partition(self) -> None:
        assert "--partition=medium" in SlurmConfig(partition="medium").resource_args()

    def test_custom_time(self) -> None:
        assert "--time=12:00:00" in SlurmConfig(time="12:00:00").resource_args()

    def test_custom_cpus(self) -> None:
        assert "--cpus-per-task=8" in SlurmConfig(cpus=8).resource_args()
