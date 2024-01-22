import os
import shutil
import logging

logging.basicConfig(level=logging.DEBUG)

trace = logging.getLogger(__name__)

def remove_file_tree(path):
    trace.debug(f"Removing file tree {path} ...")
    if os.path.exists(path):
        shutil.rmtree(path)
        trace.debug("Ok.")
    else:
        trace.debug("Not found.")

def remove_file(path):
    trace.debug(f"Removing file {path} ...")
    if os.path.exists(path):
        os.remove(path)
        trace.debug("Ok.")
    else:
        trace.debug("Not found.")

trace.info("Removing unused capd utilities...")
remove_file('src/capd_utils/capd_utils/arctan2_map.hpp')
remove_file('src/capd_utils/capd_utils/covering_factors.hpp')
remove_file('src/capd_utils/capd_utils/eigenpair.hpp')
remove_file('src/capd_utils/capd_utils/eigenproblem.hpp')
remove_file('src/capd_utils/capd_utils/eigenproblem.searcher.hpp')
remove_file('src/capd_utils/capd_utils/fixed_point_condition.hpp')
remove_file('src/capd_utils/capd_utils/ggs.hpp')
remove_file('src/capd_utils/capd_utils/grid_map.hpp')
remove_file('src/capd_utils/capd_utils/krawczyk_method.expander.hpp')
remove_file('src/capd_utils/capd_utils/lambda_map.hpp')
remove_file('src/capd_utils/capd_utils/number_string_info.cpp')
remove_file('src/capd_utils/capd_utils/number_string_info.hpp')
remove_file('src/capd_utils/capd_utils/number_string_utilities.cpp')
remove_file('src/capd_utils/capd_utils/number_string_utilities.hpp')
remove_file('src/capd_utils/capd_utils/pi.hpp')
remove_file('src/capd_utils/capd_utils/polar_coordinates.hpp')
remove_file('src/capd_utils/capd_utils/progress_logger.hpp')
remove_file('src/capd_utils/capd_utils/readable_mpreal.hpp')
remove_file('src/capd_utils/capd_utils/scalar_parser.hpp')
remove_file('src/capd_utils/capd_utils/timemap_wrapper_extended.hpp')
remove_file('src/capd_utils/capd_utils/type_adapter_map.hpp')

remove_file('src/capd_utils/capd_utils/parallel_shooting/psm.hpp')
remove_file('src/capd_utils/capd_utils/parallel_shooting/cpsm.hpp')
remove_file('src/capd_utils/capd_utils/parallel_shooting/ecpsm.hpp')
remove_file('src/capd_utils/capd_utils/parallel_shooting/epsm.hpp')
remove_file('src/capd_utils/capd_utils/parallel_shooting/epsmr.hpp')

remove_file('src/capd_utils/capd_utils/capd/fenv_rounding.hpp')

trace.info("Removing plotting related code...")
remove_file_tree('src/plot_common')
remove_file_tree('src/plot_1')
remove_file_tree('src/pcr3bp_obsolete')
remove_file_tree('src/bluesky')
remove_file_tree('src/tools/plotting')
remove_file('src/main.cpp')
remove_file('src/capd_renderable.hpp')

trace.info("Removing basic test code...")
remove_file_tree('src/pcr3bp_basic_test')

trace.info("Removing homoclinic initial test...")
remove_file('src/proof/homoclinic_orbit_origins_initial_test.cpp')

trace.info("Removing extended periodic orbit parameters test...")
remove_file('src/proof/periodic_orbit_parameters_test.extended.cpp')

trace.info("Removing collision manifold derivative check...")
remove_file('src/proof/covering_relations_test.collision_manifold_derivative_check.cpp')
remove_file('src/proof/covering_relations_test.collision_manifold_derivative_check.hpp')

trace.info("Removing covering relations setup export...")
remove_file('src/proof/covering_relations_setup_export.cpp')

trace.info("Removing print bootstrap utility...")
remove_file('src/tools/print_bootstrap.hpp')
remove_file('src/tools/print_bootstrap.cpp')


trace.info("Replacing current CMakeLists.txt ...")
remove_file('CMakeLists.txt')
shutil.copyfile('CMakeLists.proof_only.txt', 'CMakeLists.txt')

