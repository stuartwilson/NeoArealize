# The main process executable, for debugging and playing with the code.
add_executable(process process.cpp)
target_link_libraries(process ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Demonstrate Figure 5, and save data
add_executable(process_2d_karb_fig5 process_2d_karb_fig5.cpp)
target_link_libraries(process_2d_karb_fig5 ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Can't make float versions of models until libmorphologica is correctly templated/overloaded for double and float.
#add_executable(processf processf.cpp)
#target_link_libraries(processf ${MORPH_LDFLAGS})
