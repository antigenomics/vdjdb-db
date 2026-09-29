"""``python -m vdjdb``, so a worker process is the same interpreter as the one that started it.

:func:`vdjdb.annotate.junction.infer` starts its slices with ``sys.executable -m vdjdb``. Going
through the console script instead would depend on it being on ``PATH``, which it is not inside
``uv run`` subprocesses or under a SLURM step.
"""
from .cli import app

app()
