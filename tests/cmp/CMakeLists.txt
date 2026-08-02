project(testCmp LANGUAGES CXX)
viennaps_add_executable(${PROJECT_NAME} "cmp.cpp")
add_dependencies(ViennaPS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
