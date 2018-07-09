# The main process executable
add_executable(process process.cpp)
target_link_libraries(process ${MORPH_LDFLAGS})

# Can't make float versions of models until libmorphologica is correctly templated/overloaded for double and float.
#add_executable(processf processf.cpp)
#target_link_libraries(processf ${MORPH_LDFLAGS})

add_executable(testhexgrid testhexgrid.cpp HexGrid.cpp)
