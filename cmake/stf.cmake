if (TARGET stf::stf)
    return()
endif()

message(STATUS "Third-party (external): creating target 'stf::stf'")

include(CPM)
CPMAddPackage(
    NAME stf
    #GITHUB_REPOSITORY qnzhou/space-time-functions
    GIT_REPOSITORY git@github.com:qnzhou/space-time-functions.git
    GIT_TAG main
)

set_target_properties(stf PROPERTIES FOLDER third_party)
