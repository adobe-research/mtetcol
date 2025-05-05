cmake_minimum_required(VERSION 3.28)

project(mtetcol)

option(MTETCOL_BUILD_TESTS "Build tests" OFF)
option(MTETCOL_BUILD_EXAMPLE "Build example" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(strong_type)
include(unordered_dense)
include(spdlog)
include(smallvector)
include(mshio)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/mtetcol/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(mtetcol STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(mtetcol PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(mtetcol
    PUBLIC
    strong_type::strong_type
    smallvector::smallvector
    PRIVATE
    unordered_dense::unordered_dense
    spdlog::spdlog
    mshio::mshio
)

add_library(mtetcol::mtetcol ALIAS mtetcol)
target_compile_features(mtetcol PRIVATE cxx_std_20)
set_target_properties(mtetcol PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    COMPILE_WARNING_AS_ERROR ON
)

if (MTETCOL_BUILD_EXAMPLE)
    include(mtet)
    include(indirect-value-lite)
    add_executable(sweep ${CMAKE_CURRENT_LIST_DIR}/example/sweep.cpp)
    target_link_libraries(sweep
        PRIVATE
        mtetcol::mtetcol
        spdlog::spdlog
        mtet::mtet
        unordered_dense::unordered_dense
        nonstd::indirect-value-lite
    )
    target_compile_features(sweep PRIVATE cxx_std_20)
endif()

if(MTETCOL_BUILD_TESTS)
    include(CTest)
    enable_testing()
    include(catch2)

    file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
    add_executable(mtetcol_tests ${TEST_FILES})
    target_link_libraries(mtetcol_tests
        PRIVATE
        mtetcol::mtetcol
        Catch2::Catch2WithMain
        spdlog::spdlog
    )
    target_include_directories(mtetcol_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR}/src)
    target_compile_definitions(mtetcol_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
    target_compile_features(mtetcol_tests PRIVATE cxx_std_20)

    catch_discover_tests(mtetcol_tests)
endif()
