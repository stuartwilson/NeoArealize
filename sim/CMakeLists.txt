# The main process executable, for debugging and playing with the code.
add_executable(process process.cpp)
target_link_libraries(process ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Demonstrate Figure 4, and save data
add_executable(process_2d_karb_fig4 process_2d_karb_fig4.cpp)
target_link_libraries(process_2d_karb_fig4 ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Demonstrate Figure 5 (emx2 off), and save data
add_executable(process_2d_karb_fig5 process_2d_karb_fig5.cpp)
target_link_libraries(process_2d_karb_fig5 ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Demonstrate Figure 6 (pax6 off), and save data
add_executable(process_2d_karb_fig6 process_2d_karb_fig6.cpp)
target_link_libraries(process_2d_karb_fig6 ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Demonstrate Figure 8 (dual fgf8), and save data
add_executable(process_2d_karb_fig8 process_2d_karb_fig8.cpp)
target_link_libraries(process_2d_karb_fig8 ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

# Demonstrate gaussian guidance molecule system
add_executable(process_2d_karb_gaussguidance process_2d_karb_gaussguidance.cpp)
target_link_libraries(process_2d_karb_gaussguidance ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})
