# The main process executable
add_executable(process process.cpp)
target_link_libraries(process ${MORPH_LDFLAGS})
