cmake_minimum_required (VERSION 2.6)
project (NN)

# The version number.
set (NN_VERSION_MAJOR 1)
set (NN_VERSION_MINOR 0)

add_definitions("-g -Wall -Wextra")

add_library(nearestNeighbour SHARED card.c node.c rect.c rtree.c split_l.c dynamicStore.c dataBlock.c grid.c nn.c)
target_link_libraries (nearestNeighbour m)

#add_executable(test_dataBlock.x test_dataBlock.c)
#target_link_libraries(test_dataBlock.x nearestNeighbour)

#add_executable(test_grid.x test_grid.c grid.c)

#add_executable(test.x test.c)
#target_link_libraries(test.x nearestNeighbour )

add_executable(test_nn.x test_nn.cpp)
target_link_libraries(test_nn.x nearestNeighbour )

