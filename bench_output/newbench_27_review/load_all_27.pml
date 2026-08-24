# PyMOL script — load all 27 newbench_27 delivered models + crystals for review.
#
# Usage (laptop-side, after juplaunch reverse tunnel to CESGA):
#     pymol load_all_27.pml
# or from PyMOL prompt:
#     @load_all_27.pml
#
# Each PDB loads as two objects:
#     <PDB>_crystal    — deposited crystal (reference)
#     <PDB>_delivered  — FRUTON's delivery (AF-filled for 1K4Y, crystal-as-is via
#                        BL-Pose fallback for the other 26)
# Objects are grouped by PDB id so you can toggle whole entries in the sidebar.

from glob import glob
import os

for pdb in [
    "1K4Y", "2XU3", "3DF9", "3I06", "4A5S", "4B6L", "4L7G", "4QB3", "4X7Q",
    "5HJS", "5HU9", "5M7U", "5OTE", "6E1Z", "6NRH", "6Z1T", "7D5B", "8A27",
    "8ELC", "8X92", "8ZP2", "9D9I", "9GJ1", "9J9D", "9OJO", "9SI4", "9Y9N",
]:
    crystal = os.path.join("crystals",  f"{pdb}_crystal.pdb")
    delivered = os.path.join("models",  f"{pdb}_delivered.pdb")
    if os.path.isfile(crystal):
        cmd.load(crystal, f"{pdb}_crystal")
    if os.path.isfile(delivered):
        cmd.load(delivered, f"{pdb}_delivered")
    cmd.group(pdb, f"{pdb}_crystal {pdb}_delivered")

# Cosmetic defaults
cmd.hide("everything", "all")
cmd.show("cartoon", "all")
cmd.color("cyan", "*_crystal")
cmd.color("magenta", "*_delivered")
cmd.show("sticks", "resn HEM+ZN+MG+FE+CU+MN+NI+CO+CA")
cmd.show("spheres", "resn ZN+MG+FE+CU+MN+NI+CO+CA and elem Zn+Mg+Fe+Cu+Mn+Ni+Co+Ca")
cmd.bg_color("white")
cmd.zoom("all", buffer=5)

print("Loaded 27 PDB entries. Toggle a group by clicking its name in the sidebar.")
print("Cyan = crystal reference. Magenta = FRUTON delivery.")
print("For AF-filled 1K4Y, magenta ≠ cyan in the gap regions.")
print("For BL-Pose fallback (26 others), magenta == cyan (crystal shipped as-is).")
