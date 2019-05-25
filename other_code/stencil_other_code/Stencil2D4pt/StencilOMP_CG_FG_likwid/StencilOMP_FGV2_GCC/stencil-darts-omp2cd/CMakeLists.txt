add_executable(stencil.darts.exe main.output.darts.cpp stencil_omp.output.darts.cpp timer.output.darts.cpp ../../darts-lib/utils.cpp)
target_link_libraries(stencil.darts.exe darts)

include_directories(../../darts-lib)
