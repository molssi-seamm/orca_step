"""An ORCA task asks for the memory its %maxcore needs, so a run sent to a queue
does not get the queue's default (OOM-killed on Owl, ChemAI 5469, 2026-10-30)."""

from orca_step.orca_base import task_resources

MB = 1_000_000


def test_memory_follows_the_last_maxcore():
    text = (
        "! CCSD(T)-F12/RI cc-pVTZ-F12\n%pal nprocs 4 end\n"
        "%maxcore 60000\n%maxcore 24000\n"
    )
    r = task_resources(4, text)
    assert r.ntasks == 4
    # 24000 MB is 0.8 of the memory per rank
    assert r.mem_per_cpu == int(24000 * MB / 0.8)


def test_the_steps_own_maxcore():
    r = task_resources(
        2, "! B3LYP def2-SVP\n%pal nprocs 2 end\n%maxcore 3200\n* xyz 0 1\n"
    )
    assert r.mem_per_cpu == int(3200 * MB / 0.8)


def test_no_maxcore_leaves_the_default():
    assert task_resources(1, "! HF STO-3G\n").mem_per_cpu is None
