import shutil

# remove plotting related code
shutil.rmtree('src/plot_common', ignore_errors=True)
shutil.rmtree('src/plot_1', ignore_errors=True)
shutil.rmtree('src/pcr3bp_obsolete', ignore_errors=True)
shutil.rmtree('src/bluesky', ignore_errors=True)
shutil.rmtree('src/tools/plotting', ignore_errors=True)
shutil.rmtree('src/main.cpp', ignore_errors=True)
shutil.rmtree('src/capd_renderable.hpp', ignore_errors=True)

# remove basic test code
shutil.rmtree('src/pcr3bp_basic_test', ignore_errors=True)

