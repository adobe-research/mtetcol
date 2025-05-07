if (TARGET mtet::mtet)
    return()
endif()

message(STATUS "Third-party (external): creating target 'mtet::mtet'")

include(CPM)
CPMAddPackage(
  NAME mtet
  GITHUB_REPOSITORY qnzhou/mtet
  GIT_TAG 62bd4bed56d8685835aee5019d4463918bf4cbe0
)
