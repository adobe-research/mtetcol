cmake_minimum_required(VERSION 3.14)

project(mtetcol)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(strong_type)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/mtetcol/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(mtetcol STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(mtetcol PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(mtetcol
    PUBLIC
    strong_type::strong_type
)
add_library(mtetcol::mtetcol ALIAS mtetcol)
target_compile_features(mtetcol PRIVATE cxx_std_20)
set_target_properties(mtetcol PROPERTIES POSITION_INDEPENDENT_CODE ON)
