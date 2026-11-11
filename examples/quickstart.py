"""BaryCon in 30 lines: compare a model with every gas fraction in r500c."""
import matplotlib.pyplot as plt
import numpy as np

import barycon as bc

col = bc.load()
print(col.summary())

# 1) select, then homogenise explicitly
gas = col.filter(quantity="f_gas", delta="500c", region="sphere").homogenize("Planck18")
for e in gas:
    print(e.id, "|", "; ".join(e.history))

# 2) optional: correct only hydrostatic masses for a bias b=0.2
gas_b = bc.Collection(e.apply_hse_bias(0.2) if e.mass_proxy == "hydrostatic" else e for e in gas)

# 3) overplot your model (here a toy power law; e.g. baccoemu.get_baryon_fractions)
ax = gas_b.plot("f_gas", universal=0.157)
M = np.geomspace(1e13, 2e15, 100)                     # Msun/h (Planck18), as on the x-axis
ax.plot(M, 0.12 * (M / 3e14) ** 0.25, "k-", lw=2, label="my model")
ax.legend(fontsize=7, ncol=2, frameon=False)
plt.savefig("quickstart.png", dpi=150)

# 4) data for a fit, plus the references to cite
tab = gas_b.filter(data_type="points").to_table("f_gas", h_units=True)
print(tab[:5])
print(gas_b.bibtex()[:500], "...")
