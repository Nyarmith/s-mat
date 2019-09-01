file(GLOB test_srcs "*.cc")

foreach(test ${test_srcs})
    get_filename_component(exe_name ${test} NAME_WE)
    add_executable(${exe_name} ${test})
    target_link_libraries(${exe_name} s-mat)
    add_test(NAME ${test} COMMAND ${exe_name})
endforeach()
