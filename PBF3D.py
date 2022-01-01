import taichi as ti
import numpy as np
from Fluid import *

def main():
    screen_resolution = 400, 400
    max_frame = 50000

    gui = ti.GUI("pbf")
    pbfFluid = PBFFluidSim3D(gui=gui)
    pbfFluid.initialize()

    for i in range(50000):
        pbfFluid.run_pbf(i+1)
        pbfFluid.move_board(i+1)
        pbfFluid.gui_show(i+1)


if __name__ == "__main__":
    main()