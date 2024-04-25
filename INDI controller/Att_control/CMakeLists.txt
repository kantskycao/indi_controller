px4_add_library(Att_control
        Att_control.cpp
        Att_control.hpp
)
target_compile_options(Att_control PRIVATE ${MAX_CUSTOM_OPT_LEVEL})
target_include_directories(Att_control PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


