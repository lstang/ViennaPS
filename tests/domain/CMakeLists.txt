project(domain LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaPS)
target_compile_definitions(${PROJECT_NAME} PRIVATE VIENNAPS_HAS_SUNDIALS)

add_dependencies(ViennaPS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
