cmake_minimum_required(VERSION 3.20)
project(ANCF_model_c__)

set(CMAKE_CXX_STANDARD 17)
set(EIGEN_DIR "C:\\eigen-3.4.0")
include_directories (${EIGEN_DIR})


add_executable(ANCF_model_c__ main.cpp beam.cpp beam.h utils.cpp utils.h force_engine.cpp force_engine.h particle.cpp particle.h system_engine.cpp system_engine.h time_march.cpp time_march.h beam_forces.cpp beam_forces.h fluid_field.cpp fluid_field.h collision_engine.cpp collision_engine.h oriented_bounding_box.cpp oriented_bounding_box.h bounding_sphere.cpp bounding_sphere.h external_load_field.h external_load_point.h)
#target_link_libraries (ANCF_model_c__ stdc++fs)
target_link_libraries (ANCF_model_c__ ${EIGEN_DIR})