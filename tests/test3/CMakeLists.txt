project(test3 LANGUAGES NONE)

# find_program(PWSH NAMES pwsh powershell)

# fixture
add_test(NAME test.tut1_ex3.unzip
   COMMAND cmake -E tar x ${CMAKE_CURRENT_SOURCE_DIR}/ex3.zip
)

set_tests_properties(test.tut1_ex3.unzip PROPERTIES FIXTURES_SETUP UNZIP)

# tutorial 1 ex3
add_test(NAME test.tut1_ex3
    COMMAND Fastmech ex3/Case1.cgn
)

set_tests_properties(test.tut1_ex3 PROPERTIES FIXTURES_REQUIRED UNZIP)
