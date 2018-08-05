# The main process executable
add_executable(process process.cpp)
target_link_libraries(process ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Can't make float versions of models until libmorphologica is correctly templated/overloaded for double and float.
#add_executable(processf processf.cpp)
#target_link_libraries(processf ${MORPH_LDFLAGS})
