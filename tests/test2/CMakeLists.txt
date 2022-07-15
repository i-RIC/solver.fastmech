project(test2 LANGUAGES NONE)

# find_program(PWSH NAMES pwsh powershell)

# fixture
add_test(NAME test.tut1_ex3_short.unzip
   COMMAND cmake -E tar x ${CMAKE_CURRENT_SOURCE_DIR}/ex3-short.zip
)

set_tests_properties(test.tut1_ex3_short.unzip PROPERTIES FIXTURES_SETUP UNZIP)

# tutorial 1 ex3 (1000 iterations)
add_test(NAME test.tut1_ex3_short
    COMMAND Fastmech ex3-short/Case1.cgn
)

set_tests_properties(test.tut1_ex3_short PROPERTIES FIXTURES_REQUIRED UNZIP)