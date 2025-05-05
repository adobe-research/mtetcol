if (TARGET mtet::mtet)
    return()
endif()

message(STATUS "Third-party (external): creating target 'mtet::mtet'")

include(CPM)
CPMAddPackage(
  NAME mtet
  GITHUB_REPOSITORY qnzhou/mtet
  GIT_TAG e964d02b51ba62ff7b432e946722de5e80aa9f05
)
