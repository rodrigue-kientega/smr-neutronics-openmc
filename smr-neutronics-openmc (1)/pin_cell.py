"""
SMR fuel pin-cell criticality model (OpenMC)
============================================

A representative PWR/SMR fuel pin cell, built in the style of the NuScale-type
assembly I modelled during my research internship (supervised by Prof. Yahya
Tayalati). It sets up a UO2 pin with helium gap, zirconium cladding and borated
water, and runs a k-eigenvalue (criticality) calculation.

Geometry and materials are a *representative* pin cell — the kind of model one
can share openly — not proprietary plant data. The full assembly study used the
same workflow scaled to a lattice.

Requirements: OpenMC (https://docs.openmc.org) and a cross-section data library
(e.g. ENDF/B-VII.1 or VIII.0) with the OPENMC_CROSS_SECTIONS environment
variable set. Run with:  python pin_cell.py

Author: Wendlasida Rodrigue Kientega
"""
import openmc

# ----------------------------------------------------------------------
# 1. Materials
# ----------------------------------------------------------------------
uo2 = openmc.Material(name="UO2 fuel 4.05 wt% enriched")
uo2.add_element("U", 1.0, enrichment=4.05)   # weight-percent U-235
uo2.add_element("O", 2.0)
uo2.set_density("g/cm3", 10.5)

helium = openmc.Material(name="Helium gap")
helium.add_element("He", 1.0)
helium.set_density("g/cm3", 0.00178)

zirc = openmc.Material(name="Zircaloy cladding")
zirc.add_element("Zr", 1.0)
zirc.set_density("g/cm3", 6.55)

# Borated water, 1235 ppm natural boron (by mass of water)
water = openmc.Material(name="Borated water 1235 ppm")
water.add_element("H", 2.0)
water.add_element("O", 1.0)
water.add_element("B", 1235e-6, "wo")        # 1235 ppm boron by weight
water.set_density("g/cm3", 0.72)             # hot operating density
water.add_s_alpha_beta("c_H_in_H2O")          # thermal scattering in water

materials = openmc.Materials([uo2, helium, zirc, water])
materials.export_to_xml()

# ----------------------------------------------------------------------
# 2. Geometry  (standard PWR pin dimensions, square pitch)
# ----------------------------------------------------------------------
r_fuel = openmc.ZCylinder(r=0.4096)   # fuel pellet radius  (cm)
r_gap = openmc.ZCylinder(r=0.4178)    # outer gap
r_clad = openmc.ZCylinder(r=0.4750)   # outer clad
pitch = 1.26                          # cm

# Reflective box -> infinite lattice (computes k_inf of the pin cell)
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type="reflective")

fuel_cell = openmc.Cell(fill=uo2, region=-r_fuel)
gap_cell = openmc.Cell(fill=helium, region=+r_fuel & -r_gap)
clad_cell = openmc.Cell(fill=zirc, region=+r_gap & -r_clad)
water_cell = openmc.Cell(fill=water, region=+r_clad & -box)

root = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, water_cell])
openmc.Geometry(root).export_to_xml()

# ----------------------------------------------------------------------
# 3. Settings  (eigenvalue / criticality)
# ----------------------------------------------------------------------
settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.batches = 150
settings.inactive = 30
settings.particles = 20000            # 2.5M+ active histories in total

# start neutrons in the fuel
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0.0, 0.0, 0.0)),
    energy=openmc.stats.Watt())

settings.export_to_xml()

# ----------------------------------------------------------------------
# 4. Thermal-flux tally (0 - 0.625 eV), as in the assembly study
# ----------------------------------------------------------------------
tally = openmc.Tally(name="thermal flux")
e_filter = openmc.EnergyFilter([0.0, 0.625])     # eV
tally.filters = [e_filter]
tally.scores = ["flux"]
openmc.Tallies([tally]).export_to_xml()

if __name__ == "__main__":
    # Requires OpenMC + cross-section data to be installed.
    openmc.run()
    # After the run, the final k-effective is printed by OpenMC and stored in
    # the statepoint file (statepoint.150.h5). Expected: k_inf ~ 1.2-1.4 for a
    # fresh, lightly-borated pin cell at these conditions.
