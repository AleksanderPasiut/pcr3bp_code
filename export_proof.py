import os
import shutil

def remove_file(path):
    if os.path.exists(path):
        os.remove(path)

# remove plotting related code
shutil.rmtree('src/plot_common', ignore_errors=True)
shutil.rmtree('src/plot_1', ignore_errors=True)
shutil.rmtree('src/pcr3bp_obsolete', ignore_errors=True)
shutil.rmtree('src/bluesky', ignore_errors=True)
shutil.rmtree('src/tools/plotting', ignore_errors=True)
remove_file('src/main.cpp')
remove_file('src/capd_renderable.hpp')

# remove basic test code
shutil.rmtree('src/pcr3bp_basic_test', ignore_errors=True)

# remove homoclinic initial test
remove_file('src/proof/homoclinic_orbit_origins_initial_test.cpp')

# replace current CMakeLists.txt
remove_file('CMakeLists.txt')
shutil.copyfile('CMakeLists.proof_only.txt', 'CMakeLists.txt')
