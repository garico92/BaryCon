# Installation

```bash
pip install git+https://github.com/garico92/barycon   # until the PyPI release
# or, from a clone
pip install -e .            # editable
pip install -e ".[dev]"     # + pytest, black, flake8, build, twine
pip install -e ".[colossus]"  # optional: colossus concentration–mass relations
```

Check:

```python
import barycon
print(barycon.__version__, barycon.__data_version__)
from barycon import BaryonFractions        # or: from baryon_fractions import BaryonFractions
fig, ax = BaryonFractions().plot_gas_fractions_M()
```

Requirements: Python ≥ 3.9, numpy, scipy, astropy, matplotlib, pyyaml, jsonschema. Data
files are installed inside the package. Uninstall with `pip uninstall barycon`.
